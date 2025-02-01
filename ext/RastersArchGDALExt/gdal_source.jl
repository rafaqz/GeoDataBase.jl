const GDAL_LOCUS = Start()

const GDAL_DIM_ORDER = (X(), Y(), Band())

# drivers supporting the gdal Create() method to directly write to disk
const GDAL_DRIVERS_SUPPORTING_CREATE = ("GTiff", "HDF4", "KEA", "netCDF", "PCIDSK", "Zarr", "MEM"#=...=#)

# order is equal to https://gdal.org/user/virtual_file_systems.html
const GDAL_VIRTUAL_FILESYSTEMS = "/vsi" .* (
    "zip",
    "tar",
    "gzip",
    "7z",
    "rar",
    "curl",
    "curl_streaming",
    "s3",
    "s3_streaming",
    "gs",
    "gs_streaming",
    "az",
    "az_streaming",
    "adls",
    "oss",
    "oss_streaming",
    "swift",
    "swift_streaming",
    "hdfs",
    "webhdfs",
    "stdin",
    "stdout",
    "mem",
    "subfile",
    "sparse",
    )


# TODO more cases of return values here, like wrapped disk arrays
RA.sourcetrait(A::AG.RasterDataset) = GDALsource()
RA.sourcetrait(A::AG.Dataset) = GDALsource()

RA.cleanreturn(A::AG.RasterDataset) = Array(A)
RA.haslayers(::GDALsource) = false

function Base.write(filename::AbstractString, ::GDALsource, A::AbstractRasterStack; kw...)
    ext = splitext(filename)[2]
    throw(ArgumentError("Cant write a RasterStack to $ext with gdal")) 
end
function Base.write(filename::AbstractString, ::GDALsource, A::AbstractRaster{T};
    force=false, 
    verbose=true, 
    write=true,
    missingval=nokw,
    scale=nokw,
    offset=nokw,
    coerce=nokw,
    eltype=Missings.nonmissingtype(T),
    f=identity,
    kw...
) where T
    RA.check_can_write(filename, force)
    write = f === identity ? write : true
    A1 = _maybe_permute_to_gdal(A)

    # Missing values
    missingval_pair = RA._write_missingval_pair(A, missingval; eltype, verbose)

    _create_with_driver(filename, dims(A1), eltype; 
        missingval=missingval_pair[1], _block_template=A1, scale, offset, verbose, kw...
    ) do dataset
        if write
            mod = RA._mod(eltype, missingval_pair, scale, offset, coerce)
            open(A1; write=true) do O
                R = RA._maybe_modify(AG.RasterDataset(dataset), mod)
                R .= parent(O)
                if hasdim(A, Band())
                    f(R)
                else
                    f(view(R, :, :, 1))
                end
            end
        end
    end
    return filename
end

function RA._open(f, ::GDALsource, filename::AbstractString; 
    write=false, 
    mod=RA.NoMod(), 
    kw...
)
    # Check the file actually exists because the GDAL error is unhelpful
    if !isfile(filename)
        # Allow gdal virtual file systems
        # the respective string is prepended to the data source,
        # e.g. /vsicurl/https://...
        if !(length(filename) >= 8 && any(startswith.(filename, GDAL_VIRTUAL_FILESYSTEMS)))
            # Check if the filename is a url
            if RA._isurl(filename)
                filename = "/vsicurl/" * filename
            else
                # Throw our own error that the file does not exist
                RA._filenotfound_error(filename)
            end
        end
    end
    flags = write ? AG.OF_UPDATE : AG.OF_READONLY
    return AG.readraster(filename; flags) do A
        RA.cleanreturn(f(RA._maybe_modify(A, mod))) 
    end
end
RA._open(f, ::GDALsource, A::AG.RasterDataset; mod=RA.NoMod(), kw...) =
    RA.cleanreturn(f(RA._maybe_modify(A, mod)))


# DimensionalData methods for ArchGDAL types ###############################

# These methods are type piracy on DimensionalData/ArchGDAL and may have to move some day

# We allow passing in crs and mappedcrs manually
function RA._dims(raster::AG.RasterDataset, crs=nokw, mappedcrs=nokw)
    gt_dims = try
        AG.getgeotransform(raster)
    catch
        GDAL_EMPTY_TRANSFORM
    end
    gt = gt_dims
    xsize, ysize = size(raster)
    nbands = AG.nraster(raster)
    bandnames = _bandnames(raster, nbands)
    band = if all(==(""), bandnames)
        Band(Categorical(1:nbands; order=ForwardOrdered()))
    else
        Band(Categorical(bandnames; order=Unordered()))
    end

    crs = crs isa NoKW ? Rasters.crs(raster) : crs
    mappedcrs = mappedcrs isa NoKW ? nothing : mappedcrs
    xy_metadata = metadata(raster)

    # Output Sampled index dims when the transformation is lat/lon aligned,
    # otherwise use Transformed index, with an affine map.
    if _isaligned(gt)
        # Get step sizes
        xstep = gt[GDAL_WE_RES]
        ystep = gt[GDAL_NS_RES] # Usually a negative number
        # Get min, max and sampling depending on AREA_OR_POINT
        if _gdalmetadata(raster.ds, "AREA_OR_POINT") == "Point"
            sampling = Points()
            xmin, xmax = gt[GDAL_TOPLEFT_X], gt[GDAL_TOPLEFT_X] + xstep * (xsize - 1)
            ymax = gt[GDAL_TOPLEFT_Y]
            ymin = gt[GDAL_TOPLEFT_Y] + ystep * (ysize - 1)
        else
            # GeoTiff uses the "pixelCorner" convention
            sampling = Intervals(GDAL_LOCUS)
            xmin, xmax = if xstep > 0
                gt[GDAL_TOPLEFT_X], gt[GDAL_TOPLEFT_X] + xstep * (xsize - 1)
            else
                gt[GDAL_TOPLEFT_X] + xstep, gt[GDAL_TOPLEFT_X] + xstep * xsize
            end
            ymax, ymin = if ystep > 0
                gt[GDAL_TOPLEFT_Y], gt[GDAL_TOPLEFT_Y] + ystep * (ysize - 1)
            else
                gt[GDAL_TOPLEFT_Y] + ystep, gt[GDAL_TOPLEFT_Y] + ystep * ysize
            end
        end

        # Define order
        xorder = xstep > 0 ? ForwardOrdered() : ReverseOrdered()
        yorder = ystep > 0 ? ForwardOrdered() : ReverseOrdered()
        # Create lookup index. LinRange is easiest always the right size after fp error
        xindex = range(; start=xmin, stop=xmax, length=xsize)
        yindex = range(; start=ymax, stop=ymin, length=ysize)

        # Define `Projected` lookups fo X and Y dimensions
        xlookup = Projected(xindex;
            order=xorder,
            span=Regular(step(xindex)),
            sampling=sampling,
            metadata=xy_metadata,
            crs=crs,
            mappedcrs=mappedcrs,
        )
        ylookup = Projected(yindex;
            order=yorder,
            sampling=sampling,
            # Use the range step as is will be different to ystep due to float error
            span=Regular(step(yindex)),
            metadata=xy_metadata,
            crs=crs,
            mappedcrs=mappedcrs,
        )
        x = X(xlookup)
        y = Y(ylookup)

        DD.format((x, y, band), map(Base.OneTo, (xsize, ysize, nbands)))
    else
        affinemap = RA.geotransform2affine(gt)
        x = X(RA.AffineProjected(affinemap; crs, mappedcrs, metadata=xy_metadata, dim=X(), paired_lookup=Base.OneTo(ysize)))
        y = Y(RA.AffineProjected(affinemap; crs, mappedcrs, metadata=xy_metadata, dim=Y(), paired_lookup=Base.OneTo(xsize)))

        DD.format((x, y, band), map(Base.OneTo, (xsize, ysize, nbands)))
    end
end

# TODO make metadata optional, its slow to get
function RA._metadata(rds::AG.RasterDataset, args...)
    band = AG.getband(rds.ds, 1)
    metadata = RA._metadatadict(GDALsource())
    # color = AG.getname(AG.getcolorinterp(band))
    scale = AG.getscale(band)
    offset = AG.getoffset(band)
    # norvw = AG.noverview(band)
    units = AG.getunittype(band)
    filepath = _getfilepath(rds)
    # Set metadata if they are not default values
    if scale != oneunit(scale) 
        metadata["scale"] = scale 
    end
    if offset != zero(offset) 
        metadata["offset"] = offset
    end
    if units != ""
        metadata["units"] = units
    end
    if !isnothing(filepath)
        metadata["filepath"] = filepath
    end
    return metadata
end

# Rasters methods for ArchGDAL types ##############################

# Create a Raster from a dataset
RA.Raster(ds::AG.Dataset; kw...) = Raster(AG.RasterDataset(ds); kw...)
function RA.Raster(ds::AG.RasterDataset; lazy=false, kw...) 
    filepath = if lazy
        fp = _getfilepath(ds)
        isnothing(fp) ? "/vsimem" : fp
    else
        ""
    end
    Raster(ds, filepath; kw...)
end

RA.missingval(ds::AG.Dataset, args...) = RA.missingval(AG.RasterDataset(ds))
function RA.missingval(rds::AG.RasterDataset, args...)
    # All bands have the same missingval in GDAL
    band = AG.getband(rds.ds, 1)
    # GDAL will set this
    hasnodataval = Ref(Cint(0))
    # Int64 and UInt64 need special casing in GDAL
    nodataval = if eltype(rds) == Int64
        AG.GDAL.gdalgetrasternodatavalueasint64(band, hasnodataval)
    elseif eltype(rds) == UInt64
        AG.GDAL.gdalgetrasternodatavalueasuint64(band, hasnodataval)
    else
        AG.GDAL.gdalgetrasternodatavalue(band, hasnodataval)
    end
    # If there was a missingval, clean it and use it.
    # Otherwise set missingval to `nothing`, meaning there isnt one
    if Bool(hasnodataval[])
        return _missingval_from_gdal(eltype(band), nodataval)
    else
        return nothing
    end
end

# GDAL always returns well known text crs
function RA.crs(raster::AG.RasterDataset, args...)
    WellKnownText(GeoFormatTypes.CRS(), string(AG.getproj(raster.ds)))
end

# ArchGDAL methods for Rasters types ####################################

# Extend ArchGDAL RasterDataset and Dataset to accept AbstractRaster as input
function AG.Dataset(f::Function, A::AbstractRaster; kw...)
    AG.RasterDataset(A; kw...) do rds
        f(rds.ds)
    end
end
function AG.RasterDataset(f::Function, A::AbstractRaster; 
    filename="", 
    scale=nokw,
    offset=nokw,
    coerce=nokw,
    verbose=false,
    eltype=Missings.nonmissingtype(eltype(A)),
    missingval=Rasters.missingval(A),
    kw...
)
    return open(_maybe_permute_to_gdal(A)) do O 
        _create_with_driver(filename, dims(A), eltype; 
            _block_template=A, missingval, scale, offset, verbose, kw...
        ) do dataset
            rds = AG.RasterDataset(dataset)
            mv = RA.missingval(rds) => RA.missingval(O)
            mod = RA._mod(eltype, mv, scale, offset, coerce)
            RA._maybe_modify(rds, mod) .= parent(O)
            f(rds)
        end
    end
end

# Utils ########################################################################

# Sometimes GDAL stores the `missingval` in the wrong type, so fix it.
_missingval_from_gdal(T::Type{<:AbstractFloat}, x::Real) = convert(T, x)
function _missingval_from_gdal(T::Type{<:Integer}, x::AbstractFloat; verbose=true)
    if trunc(x) === x && x >= typemin(T) && x <= typemax(T)
        convert(T, x)
    else
        verbose && @warn "Missing value $x can't be converted to array eltype $T. `missingval` set to `nothing`"
        nothing
    end
end
function _missingval_from_gdal(T::Type{<:Integer}, x::Integer; verbose=true)
    if x >= typemin(T) && x <= typemax(T)
        convert(T, x)
    else
        verbose && @warn "Missing value $x can't be converted to array eltype $T. `missingval` set to `nothing`"
        nothing
    end
end
_missingval_from_gdal(T, x) = x

# Make sure driver is sensible
function _check_driver(::Nothing, driver) 
    isnokwornothing(driver) || isempty(driver) ? "MEM" : driver
end
function _check_driver(filename::AbstractString, driver)
    if isnokwornothing(driver) || isempty(driver)
        if isempty(filename)
            driver = "MEM"
        else
            driver = AG.extensiondriver(filename)
            if driver == "COG"
                driver = "GTiff"
            end
        end
    end
    return driver
end

# Handle creating a dataset with any driver,
# applying the function `f` to the created dataset
function _create_with_driver(f, filename, dims::Tuple, T;
    verbose=true,
    missingval=nokw,
    options=nokw,
    driver=nokw, 
    chunks=nokw,
    scale=nokw,
    offset=nokw,
    _block_template=nothing, 
    kw...
)
    # Allow but discourage south-up
    verbose && _maybe_warn_south_up(dims, verbose, "Creating a South-up raster. You may wish to reverse the `Y` dimension to use conventional North-up")
    options = isnokwornothing(options) ? Dict{String,String}() : options

    # Pairs should not get this far
    @assert !(missingval isa Pair)
    # If missingval is missing, generate one GDAL can write
    missingval = ismissing(missingval) ? RA._writeable_missing(T; verbose) : missingval
    # Make sure dimensions are valid for GDAL
    _gdal_validate(dims)

    # Move x and y locus to start
    x, y = map(DD.dims(dims, (XDim, YDim))) do d
        maybeshiftlocus(Start(), RA.nolookup_to_sampled(d))
    end

    # Band handling - a 2d raster has 1 band
    newdims = hasdim(dims, Band()) ? (x, y, DD.dims(dims, Band)) : (x, y)
    nbands = hasdim(dims, Band) ? length(DD.dims(dims, Band())) : 1

    # Driver detection
    driver = _check_driver(filename, driver)
    options_vec = _process_options(driver, options; _block_template, chunks, verbose)
    gdaldriver = driver isa String ? AG.getdriver(driver) : driver

    # Keywords and filenames for AG.create
    create_kw = (; width=length(x), height=length(y), nbands, dtype=T,)
    filename = isnothing(filename) ? "" : filename

    # Not all drivers can be directly created, some need an intermediate step
    if AG.shortname(gdaldriver) in GDAL_DRIVERS_SUPPORTING_CREATE
        AG.create(filename; driver=gdaldriver, options=options_vec, create_kw...) do dataset
            _set_dataset_properties!(dataset, newdims, missingval, scale, offset)
            f(dataset)
        end
    else
        # Create a tif and copy it to `filename`, as ArchGDAL.create
        # does not support direct creation of ASCII etc. rasters
        tif_options_vec = _process_options("GTiff", Dict{String,String}(); 
            chunks, _block_template, verbose
        )
        tif_driver = AG.getdriver("GTiff")
        tif_name = tempname() * ".tif"
        AG.create(tif_name; driver=tif_driver, options=tif_options_vec, create_kw...) do dataset
            _set_dataset_properties!(dataset, newdims, missingval, scale, offset)
            f(dataset)
            target_ds = AG.copy(dataset; filename=filename, driver=gdaldriver, options=options_vec)
            AG.destroy(target_ds)
        end
    end
end

# Makie sure gdal can actually write this raster
@noinline function _gdal_validate(dims)
    all(hasdim(dims, (XDim, YDim))) || throw(ArgumentError("`Raster` must have both an `X` and `Y` to be converted to an ArchGDAL `Dataset`"))
    if length(dims) === 3
        otherdim = otherdims(dims, (XDim, YDim))[1]
        otherdim isa Band || throw(ArgumentError("ArchGDAL can't handle $(basetypeof(thirddim)) dims - only X, Y, and Band"))
    elseif !(length(dims) in (2, 3))
        throw(ArgumentError("ArchGDAL can only accept 2 or 3 dimensional arrays"))
    end
end

# Convert a Dict of options to a Vector{String} for GDAL
function _process_options(driver::String, options::Dict; 
    chunks=nokw,
    _block_template=nothing,
    verbose=true,
)
    options_str = Dict(string(k)=>string(v) for (k,v) in options)
    # Get the GDAL driver object
    gdaldriver = AG.getdriver(driver)

    # set default compression
    if driver != "MEM" && !("COMPRESS" in keys(options_str)) && AG.validate(gdaldriver, ["COMPRESS=ZSTD"])
        options_str["COMPRESS"] = "ZSTD"
    end

    # the goal is to set write block sizes that correspond to eventually blocked reads
    # creation options are driver dependent

    chunk_pattern = RA._chunks_to_tuple(_block_template, (X(), Y(), Band()), chunks) 
    if !isnothing(chunk_pattern)
        xchunksize, ychunksize = chunk_pattern
        
        block_x, block_y = string.((xchunksize, ychunksize))

        if driver == "GTiff"
            # GDAL default is line-by-line compression without tiling.
            # Here, tiling is enabled if the source chunk size is viable for GTiff,
            # i.e. when the chunk size is divisible by 16.
            if (xchunksize % 16 == 0) && (ychunksize % 16 == 0)
                options_str["TILED"] = "YES"
            else
                xchunksize == 1 || (verbose && @warn "X and Y chunk size do not match. Columns are used and X size $xchunksize is ignored")
            end
            # don't overwrite user specified values
            if !("BLOCKXSIZE" in keys(options_str))
                options_str["BLOCKXSIZE"] = block_x
            end
            if !("BLOCKYSIZE" in keys(options_str))
                options_str["BLOCKYSIZE"] = block_y
            end
        elseif driver == "COG"
            if !("BLOCKSIZE" in keys(options_str))
                if xchunksize == ychunksize 
                    options_str["BLOCKSIZE"] = block_x
                else
                    verbose && @warn "Writing COG X and Y chunks do not match: $block_x, $block_y. Default of 512, 512 used."
                end
            end
        end
    end
    # if the input is unchunked we just use the driver defaults
    options_vec = ["$(uppercase(k))=$(uppercase(v))" for (k,v) in options_str]

    invalid_options = String[]
    for option in options_vec
        if !AG.validate(gdaldriver, [option])
            push!(invalid_options, option)
        end
    end

    if length(invalid_options) > 0
        throw(ArgumentError("Invalid driver creation option(s) detected.
        Please check them carefully with the documentation at https://gdal.org/drivers/raster/index.html.
        $invalid_options"))
    end

    return options_vec
end

# Bands can have text names, but usually don't
function _bandnames(rds::AG.RasterDataset, nbands=AG.nraster(rds))
    map(1:nbands) do b
        AG.getband(rds.ds, b) do band
            AG.GDAL.gdalgetdescription(band.ptr)
        end
    end
end

function _gdalmetadata(dataset::AG.Dataset, name)
    meta = AG.metadata(dataset)
    regex = Regex("$name=(.*)")
    i = findfirst(f -> occursin(regex, f), meta)
    if i isa Nothing
        return ""
    else
        return match(regex, meta[i])[1]
    end
end

# Set the properties of an ArchGDAL Dataset to match
# the dimensions and missingval of a Raster
_set_dataset_properties!(ds::AG.Dataset, A, scale, offset) =
    _set_dataset_properties!(ds, dims(A), missingval(A), scale, offset)
function _set_dataset_properties!(dataset::AG.Dataset, dims::Tuple, missingval, scale, offset)
    # We cant write mixed Points/Intervals, so default to Intervals if mixed
    xy = DD.dims(dims, (X, Y))
    if any(x -> x isa Intervals, map(sampling, xy)) && any(x -> x isa Points, map(sampling, xy))
        dims = set(dims, X => Intervals, Y => Intervals)
    end
    # Convert the dimensions to `Projected` if they are `Converted`
    # This allows saving NetCDF to Tiff.
    x = convertlookup(Projected, DD.dims(dims, X))
    y = convertlookup(Projected, DD.dims(dims, Y))

    # Set the index loci to the start of the cell for the lat and lon dimensions.
    # NetCDF or other formats use the center of the interval, so they need conversion.
    x = DD.maybeshiftlocus(GDAL_LOCUS, x)
    y = DD.maybeshiftlocus(GDAL_LOCUS, y)

    # Set GDAL AREA_OR_POINT metadata
    area_or_point = sampling(x) isa Points ? "Point" : "Area"
    AG.GDAL.gdalsetmetadataitem(dataset, "AREA_OR_POINT", area_or_point, "")

    # Set crs if it exists, converting crs to WKT
    if !isnothing(crs(x))
        AG.setproj!(dataset, convert(String, convert(WellKnownText, crs(x))))
    end

    # Set the geotransform from the updated lookups
    gt = RA.dims2geotransform(x, y)
    AG.setgeotransform!(dataset, gt)

    if !RA.isnokwornothing(missingval)
        bands = hasdim(dims, Band) ? axes(DD.dims(dims, Band), 1) : 1
        for i in bands
            rasterband = AG.getband(dataset, i)
            (RA.isnokw(offset) || isnothing(offset)) || AG.setoffset!(rasterband, offset)
            (RA.isnokw(scale) || isnothing(scale)) || AG.setscale!(rasterband, scale)
            # Set the missing value/nodataval. This is a little complicated
            # because gdal has separate method for 64 bit integers
            if missingval isa Int64
                AG.GDAL.gdalsetrasternodatavalueasint64(rasterband, missingval)
            elseif missingval isa UInt64
                AG.GDAL.gdalsetrasternodatavalueasuint64(rasterband, missingval)
            else
                AG.GDAL.gdalsetrasternodatavalue(rasterband, missingval)
            end
        end
    end

    # Write band labels if they are not Integers.
    if hasdim(dims, Band)
        bandlookup = DD.lookup(DD.dims(dims, Band))
        if !(eltype(bandlookup) <: Integer)
            for i in eachindex(bandlookup)
                AG.getband(dataset, i) do band
                    AG.GDAL.gdalsetdescription(band.ptr, string(bandlookup[i]))
                end
            end
        end
    end

    return dataset
end

# Detect the driver string from the extension
_extensiondriver(filename::Nothing) = "MEM"
function _extensiondriver(filename::AbstractString)
    # TODO move this check to ArchGDAL
    if filename == "/vsimem/tmp"
        "MEM"
    elseif splitext(filename)[2] == ".tif"
        # Force GTiff as the default for .tif because COG cannot do `create` yet
        "GTiff"
    else
        AG.extensiondriver(filename)
    end
end

# Permute dims unless they match the normal GDAL dimension order
_maybe_permute_to_gdal(A) = _maybe_permute_to_gdal(A, DD.dims(A, GDAL_DIM_ORDER))
_maybe_permute_to_gdal(A, dims::Tuple) = A
_maybe_permute_to_gdal(A, dims::Tuple{<:XDim,<:YDim,<:Band}) = permutedims(A, dims)
_maybe_permute_to_gdal(A, dims::Tuple{<:XDim,<:YDim}) = permutedims(A, dims)

_maybe_restore_from_gdal(A, dims::Tuple) = _maybe_reorder(permutedims(A, dims), dims)
_maybe_restore_from_gdal(A, dims::Union{Tuple{<:XDim,<:YDim,<:Band},Tuple{<:XDim,<:YDim}}) =
    _maybe_reorder(A, dims)

function _maybe_reorder(A, dims::Tuple)
    if all(map(l -> l isa AbstractSampled, lookup(dims, (XDim, YDim)))) &&
        all(map(l -> l isa AbstractSampled, lookup(A, (XDim, YDim))))
        reorder(A, dims)
    else
        A
    end
end

function _maybe_warn_south_up(A, verbose, msg)
    if hasdim(A, Y())
        verbose && lookup(A, Y()) isa AbstractSampled && order(A, Y()) isa ForwardOrdered && @warn msg
    end
    return nothing
end

#= Geotranforms ########################################################################

See https://lists.osgeo.org/pipermail/gdal-dev/2011-July/029449.html

"In the particular, but common, case of a “north up” image without any rotation or
shearing, the georeferencing transform takes the following form" :
adfGeoTransform[0] /* top left x */
adfGeoTransform[1] /* w-e pixel resolution */
adfGeoTransform[2] /* 0 */
adfGeoTransform[3] /* top left y */
adfGeoTransform[4] /* 0 */
adfGeoTransform[5] /* n-s pixel resolution (negative value) */
=#

# These function are defined in ext/RastersCoordinateTransformationsExt.jl
const USING_COORDINATETRANSFORMATIONS_MESSAGE =
    "Run `using CoordinateTransformations` to load affine transformed rasters"

function RA.dims2geotransform(x::XDim, y::YDim)
    sampling(x) == sampling(y) || throw(ArgumentError("Sampling of x and y must match to generate a geotransform"))
    (order(x) isa Ordered && order(y) isa Ordered) || throw(ArgumentError("GDAL can only write ordered lookups"))
    (span(x) isa Regular && span(y) isa Regular) || throw(ArgumentError("GDAL can only write regular lookups"))
    gt = zeros(6)
    gt[GDAL_ROT1] = zero(eltype(gt))
    gt[GDAL_ROT2] = zero(eltype(gt))
    gt[GDAL_WE_RES] = step(x)
    gt[GDAL_NS_RES] = step(y)
    if sampling(x) isa Points
        gt[GDAL_TOPLEFT_X] = first(x)
        gt[GDAL_TOPLEFT_Y] = first(y)
    else
        sampling(x) isa Intervals{Start} || throw(ArgumentError("GDAL can only write intervals with Start locus"))
        gt[GDAL_TOPLEFT_X] = order(x) isa ReverseOrdered ? first(x) - step(x) : first(x)
        gt[GDAL_TOPLEFT_Y] = order(y) isa ReverseOrdered ? first(y) - step(y) : first(y)
    end
    return gt
end
RA.geotransform2affine(gt) = error(USING_COORDINATETRANSFORMATIONS_MESSAGE)
RA.affine2geotransform(am) = error(USING_COORDINATETRANSFORMATIONS_MESSAGE)

_isaligned(geotransform) = geotransform[GDAL_ROT1] == 0 && geotransform[GDAL_ROT2] == 0

function _getfilepath(ds)
    filelist = AG.filelist(ds)
    if length(filelist) == 0
        return nothing
    else
        return first(filelist)
    end
end

# precompilation
# function _precompile(::Type{GDALsource})
#     ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

#     for T in (Any, UInt8, UInt16, Int16, UInt32, Int32, Float32, Float64)
#         DS = AG.RasterDataset{T,AG.Dataset}
#         precompile(crs, (DS,))
#         precompile(Rasters.FileArray, (DS, String))
#         precompile(dims, (DS,))
#         precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},Nothing))
#         precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},EPSG))
#         precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},ProjString))
#         precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},WellKnownText{GeoFormatTypes.CRS}))
#         precompile(metadata, (DS, key))
#         precompile(missingval, (DS, key))
#         precompile(Raster, (DS, key))
#         precompile(Raster, (DS, String, Nothing))
#         precompile(Raster, (DS, String, Symbol))
#     end
# end

# _precompile(GRDsource)

