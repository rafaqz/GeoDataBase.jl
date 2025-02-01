
const SOURCE_WRITE_DOCSTRING = """
Other keyword arguments are passed to the `write` method for the backend.

## NetCDF keywords

- `append`: If true, the variable of the current Raster will be appended
    to `filename`, if it actually exists.
- `deflatelevel`: Compression level: `0` (default) means no compression and `9`
    means maximum compression. Each chunk will be compressed individually.
- `shuffle`: If `true`, the shuffle filter is activated which can improve the
    compression ratio.
- `checksum`: The checksum method can be `:fletcher32` or `:nochecksum`,
    the default.
- `typename`: The name of the NetCDF type required for vlen arrays
    (https://web.archive.org/save/https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-c/nc_005fdef_005fvlen.html)

## GDAL Keywords

$FORCE_KEYWORD
- `driver`: A GDAL driver name `String` or a GDAL driver retrieved via `ArchGDAL.getdriver(drivername)`.
    By default `driver` is guessed from the filename extension.

- `options::Dict{String,String}`: A dictionary containing the dataset creation options passed to the driver.
    For example: `Dict("COMPRESS" => "DEFLATE")`.

Valid `driver` names and the `options` for each can be found at:
[https://gdal.org/drivers/raster/index.html](https://gdal.org/drivers/raster/index.html)


## Source comments

### R grd/grid files

Write a `Raster` to a .grd file with a .gri header file.
Returns the base of `filename` with a `.grd` extension.

### GDAL (tiff, and everything else)

Used if you `write` a `Raster` with a `filename` extension that no other backend can write.
GDAL is the fallback, and writes a lot of file types, but is not guaranteed to work.
"""

"""
    Base.write(filename::AbstractString, A::AbstractRaster; [source], kw...)

Write an [`AbstractRaster`](@ref) to file, guessing the backend from the
file extension or using the `source` keyword.

# Keywords

$CHUNKS_KEYWORD
$FORCE_KEYWORD
$WRITE_MISSINGVAL_KEYWORD
$SOURCE_KEYWORD

$SOURCE_WRITE_DOCSTRING

Returns `filename`.
"""
function Base.write(filename::AbstractString, A::AbstractRaster;
    source=sourcetrait(filename),
    missingval=nokw,
    kw...
)
    missingval = isnokw(missingval) ? Rasters.missingval(A) : missingval
    write(filename, sourcetrait(source), A; missingval, kw...)
end
Base.write(A::AbstractRaster; kw...) = write(filename(A), A; kw...)
# Fallback
function Base.write(
    filename::AbstractString, source::Source, A::Union{AbstractRaster,AbstractRasterStack}; kw...
)
    missing_package = SOURCE2PACKAGENAME[source]
    error("Missing package extension for $source. Run `using $missing_package` before using `write` for this file extension.")
end

"""
    Base.write(filename::AbstractString, s::AbstractRasterStack; kw...)

Write any [`AbstractRasterStack`](@ref) to one or multiple files,
depending on the backend. Backend is guessed from the filename extension
or forced with the `source` keyword.

If the source can't be saved as a stack-like object, individual array layers will be saved.

## Keywords

$CHUNKS_KEYWORD
$EXT_KEYWORD
$FORCE_KEYWORD
$MISSINGVAL_KEYWORD
    For `RasterStack` this may be a `NamedTuple`, one for each layer.
$SOURCE_KEYWORD
$SUFFIX_KEYWORD
$VERBOSE_KEYWORD

$SOURCE_WRITE_DOCSTRING
"""
function Base.write(path::AbstractString, s::AbstractRasterStack;
    suffix=nothing,
    ext=nothing,
    source=sourcetrait(path, ext),
    verbose=true,
    missingval=nokw,
    kw...
)
    source = sourcetrait(source)
    missingval = _stack_missingvals(s, missingval)
    if haslayers(source)
        write(path, source, s; missingval, kw...)
    else
        # Otherwise write separate files for each layer
        if isnothing(ext)
            base, ext = splitext(path)
        else
            base, _ = splitext(path)
        end
        suffix1 = if isnothing(suffix)
            divider = Sys.iswindows() ? '\\' : '/'
            # Add an underscore to the key if there is a file name already
            spacer = last(path) == divider ? "" : "_"
            map(k -> string(spacer, k), keys(s))
        else
            suffix
        end
        if verbose
            @warn string("Cannot write complete stacks to \"", ext, "\", writing layers as individual files")
        end
        map(keys(s), suffix1, missingval) do key, suf, mv
            fn = string(base, suf, ext)
            write(fn, source, s[key]; missingval=mv, verbose, kw...)
        end |> NamedTuple{keys(s)}
    end
end

_stack_missingvals(::RasterStack{<:Any,T}, x) where T = _stack_missingvals(T, x)
function _stack_missingvals(::Type{T}, missingval::NamedTuple{K}) where {K,T<:NamedTuple{K}}
    map(_types(T), values(missingval)) do t, mv
        ismissing(mv) ? _type_missingval(t) : mv
    end |> NamedTuple{K}
end
_stack_missingvals(::Type{T}, missingval::NamedTuple{K1}) where {K1,T<:NamedTuple{K2}} where K2 =
    throw(ArgumentError("stack keys $K1 do not match misssingval keys $K2"))
_stack_missingvals(::Type{T}, missingval::Missing) where T<:NamedTuple{K} where K =
    NamedTuple{K}(map(_type_missingval, _types(T)))
_stack_missingvals(::Type{T}, missingval) where T =
    _stack_nt(T, missingval)

_stack_nt(::NamedTuple{K}, x) where K = NamedTuple{K}(map(_ -> x, K))
_stack_nt(::RasterStack{<:Any,T}, x) where T = _stack_nt(T, x)
_stack_nt(::Type{T}, x::NamedTuple{K}) where {K,T<:NamedTuple{K}} = x
_stack_nt(::Type{T}, x::NamedTuple{K1}) where {K1,T<:NamedTuple{K2}} where K2 =
    throw(ArgumentError("stack keys $K1 do not match misssingval keys $K2"))
_stack_nt(::Type{T}, x) where T<:NamedTuple{K} where K =
    NamedTuple{K}(map(_ -> x, K))

@generated function _types(::Type{<:NamedTuple{K,T}}) where {K,T}
    Tuple(T.parameters)
end

"""
    Base.write(filepath::AbstractString, s::AbstractRasterSeries; kw...)

Write any [`AbstractRasterSeries`](@ref) to multiple files,
guessing the backend from the file extension.

The lookup values of the series will be appended to the filepath (before the extension),
separated by underscores.

All keywords are passed through to these `Raster` and `RasterStack` methods.

## Keywords

$CHUNKS_KEYWORD
$EXT_KEYWORD
$FORCE_KEYWORD
$MISSINGVAL_KEYWORD
    For series with `RasterStack` child objects, this may be a `NamedTuple`, one for each layer.
$SOURCE_KEYWORD
$VERBOSE_KEYWORD
"""
function Base.write(path::AbstractString, A::AbstractRasterSeries;
    ext=nothing,
    source=sourcetrait(path, ext),
    verbose=true,
    kw...
)
    source = sourcetrait(source)
    base, path_ext = splitext(path)
    ext = isnothing(ext) ? path_ext : ext
    map(A, DimPoints(A)) do raster, p
        lookupstring = join(map(string, p), "_")
        written_paths = if raster isa RasterStack && !haslayers(source)
            stack_dir = joinpath(dirname(base), lookupstring)
            mkpath(stack_dir)
            stack_filepath = joinpath(stack_dir, basename(path))
            write(stack_filepath, raster; source, verbose, ext, kw...)
        else
            dirpath = dirname(base)
            filename = basename(base) == "" ? lookupstring : join((basename(base), lookupstring), '_')
            slice_filename = joinpath(dirpath, string(filename, ext))
            write(slice_filename, raster; source, verbose, kw...)
        end
        # Don't print warnings for every slice of the series
        verbose = false
        written_paths
    end
end
Base.write(f::Function, args...; kw...) = write(args...; f, kw...)

# Trait for source data that has stack layers
haslayers(T) = false

#  This is used in source `write` methods
check_can_write(filename::Union{Nothing,NoKW}, force) = true
function check_can_write(filename, force)
    if !check_can_write(Bool, filename, force)
        throw(ArgumentError("filename already exists at $filename. use the keyword `force=true` to write anyway"))
    end
    return true
end
check_can_write(::Type{Bool}, filename::Union{Nothing,NoKW}, force) = true
check_can_write(::Type{Bool}, filename, force) = (force || (!isfile(filename) && !isdir(filename)))
