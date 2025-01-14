const TypeNamedTuple = NamedTuple{<:Any,<:Tuple{Vararg{Type}}}

"""
    create([f!], [filename], template; kw...)

Create a new, uninitialised [`Raster`](@ref) or [`RasterStack`](@ref).

If `filename` is a `String` it will be created on disk, and opened lazily.
If it is `nothing` or not passed, an in-memory `Raster` will be created.

If type is a `Type` return value is a `Raster`. The `eltype` will usually be `T`, except
where `scale` and/or `offset` keywords are used or a `missingval` of a different type is specified, 
in which case `T` will depend on the type promotion of `scale`, `offset` and `missingval` with `T`.
If `missingval` is a `Pair` of `on_disk_missingval => user_facing_missingval`, the user facing value
will effect `T`, not the internal on-disk value.

If types is a `NamedTuple` of types, the result will be a `RasterStack`. In this case `fill` and 
`missingval` can be single values (for all layers) or `NamedTuple` with the same names to specify per-layer.

`f!` will be applied to the `Raster` or `RasterStack` while it is stil open after creation, 
to avoid opening it twice. The return value of `f!` is disregarded but modifications
to the `Raster` or the `RasterStack` layers will be written to disk or changed in memory.

## Arguments

- `filename`: a String file path, which will create a file on disk and return it as
    a lazy `Raster`, or `nothing` to create an in-memory `Raster`.
- `template`: a `Raster`, `Tuple` of `Dimension` or `Extents.Extent` to use as a template.
    If an `Extent` is used, a `size` or `res` keyword must be passed.
    If a `T` argument is not used, it is taken from the `template` eltype.
- `type`: the element type to use in the created array. A `NamedTuple` of types
    will create a `RasterStack`

## Keywords

$NAME_KEYWORD
$REFDIMS_KEYWORD
$METADATA_KEYWORD
$WRITE_MISSINGVAL_KEYWORD
- `fill`: A value to fill the array with, before `scale` and `offset` are applied. 
    If there is no `fill`, raster values may remain undefined. They may be set to 
    `missingval` on disk, but this is not guaranteed. It us often more efficient to 
    use `fill` than to fill manually after `create`.
$SOURCE_KEYWORD
- `lazy`: A `Bool` specifying if to load data lazily from disk. For `create`
    `lazy=true` is the default, as creating a disk-based file is normally associated
    with it being larger than memory.
$CHUNKS_KEYWORD
$SCALE_KEYWORD
$OFFSET_KEYWORD
$COERCE_KEYWORD
$VERBOSE_KEYWORD
$RES_KEYWORD
$SIZE_KEYWORD
$CRS_KEYWORD
$CHUNKS_KEYWORD
- `reverse_y`: often we want to write `Y` dimensions in reverse.
    When building dimensions from an `Extents.Extent` and `size` or `res` we can do this by
    using `reverse_y=true`. Using a negative value in `res` will acheive the same result.
    With a template `Raster` or a `Tuple` of `Dimension`, the existing order is used.

## Example

Here we create a UInt8 GeoTIFF and open it as a Raster, from -80 to 80
lattitude, and 0 to 120 longitude, with a resolution of 0.25 degrees.

We scale values from 0-1 over `UInt8` 0-200, and using `255`.
Values that don't convert exactly will error (we could use `coerce=trunc` to fix that).

We use `UInt8(255)` as the `missingval` on disk, but mask it with `missing` in
the loaded `Raster`.

We use standard lat/lon (EPSG:4326) as the crs, and force writing if the file exists.

```julia
using Rasters, NCDatasets, ArchGDAL, Extents, Dates
using Rasters.Lookups
rast = Rasters.create("created.tif", UInt8, Extents.Extent(X=(0, 120), Y=(-80, 80), Band=(0, 12));
    res=(X=10.0, Y=10.0, Band=1),
    # size=(X=100, Y=100, Band=12),
    name=:myraster,
    crs=EPSG(4326),
    force=true,
    fill=0x01,
    sampling=(X=Intervals(Start()), Y=Intervals(Start()), Band=Intervals(Start())),
) do A
    # While we have the newly created raster open, we can write to it
    A[X=1:10, Y=1:10] .= 0xff
end

read(rast)
```

We can also create a `RasterStack` by passing a `NamedTuple` of types:

```julia
ext = Extents.Extent(X=(0, 120), Y=(-80, 80))#, Band=(1, 3))
types = (a=UInt8, b=Int32, c=Float64)
rast = Rasters.create("created.nc", types, ext;
    # res=(X=1.0, Y=1.0, Band=1),
    size=(X=100, Y=100),
    crs=EPSG(4326),
    force=true,
    # sampling=(X=Intervals(Start()), Y=Intervals(Start()), Band=Points()),
end

RasterStack("created.nc")

╭───────────────────────────────────────────╮
│ 480×640 Raster{Union{Missing, Float64},2} │
├───────────────────────────────────────────┴───────────────────────────────────────── dims ┐
  ↓ X Projected{Float64} LinRange{Float64}(0.0, 119.75, 480) ForwardOrdered Regular Points,
  → Y Projected{Float64} LinRange{Float64}(79.75, -80.0, 640) ReverseOrdered Regular Points
├───────────────────────────────────────────────────────────────────────────────── metadata ┤
  Metadata{GDALsource} of Dict{String, Any} with 2 entries:
  "filepath" => "created.tif"
  "scale"    => 0.005
├─────────────────────────────────────────────────────────────────────────────────── raster ┤
  extent: Extent(X = (0.0, 119.75), Y = (-80.0, 79.75))
  missingval: missing
  crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199
433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
  filename: nothing
└───────────────────────────────────────────────────────────────────────────────────────────┘
```
"""
# Create with a function that will be called to fill the raster
create(f::Function, args...; kw...) = create(args...; kw..., f)
# Create from Raster or RasterStack with no filename
create(A::Union{AbstractRaster,AbstractRasterStack}; kw...) = create(nothing, A; kw...)
# Create from type and Raster or RasterStack with no filename
create(T::Union{Type,TypeNamedTuple}, A::Union{Tuple,Extents.Extent,AbstractRaster,AbstractRasterStack}; kw...) =
    create(nothing, T, A; kw...)
# Create from filename and Raster
function create(filename::Union{AbstractString,Nothing}, A::AbstractRaster{T}; 
    missingval=missingval(A), # Only take missingval here when types are not specified
    kw...
) where T
    create(filename, T, A; missingval, kw...)
end
# Create from filename and RasterStack
function create(filename::Union{AbstractString,Nothing}, st::AbstractRasterStack; 
    missingval=missingval(st), # Only take missingval here when types are not specified
    kw...
)
    create(filename, map(eltype, layers(st)), st; missingval, kw...)
end
# Create Raster from filename, type and a Raster, using its
# parent as the parent type so e.g. CuArray will propagate
create(filename::Union{AbstractString,Nothing}, T::Union{Type,TypeNamedTuple}, A::AbstractRaster; kw...) =
    create(filename, T, dims(A); parent=parent(A), kw...)
# Create RasterStack from filename, NamedTuple type and dims
function create(filename::Union{AbstractString,Nothing}, T::NamedTuple{K1}, st::AbstractRasterStack{K2};
    metadata=metadata(st),
    layerdims=nokw,
    layermetadata=nokw,
    kw...
) where {K1,K2}
    if all(map(in(K2), K1))
        layerdims = isnokw(layerdims) ? DD.layerdims(st)[K1] : layerdims
        layermetadata = isnokw(layermetadata) ? DD.layermetadata(st)[K1] : layermetadata
    end
    return create(filename, T, dims(st);
        parent=first(parent(st)), metadata, missingval, layerdims, layermetadata, kw...
    )
end
# Create from filename, type and dims
function create(filename::AbstractString, T::Union{Type,NamedTuple}, dims::Tuple;
    lazy=true,
    parent=nokw,
    suffix=nokw,
    source::Source=_sourcetrait(filename),
    missingval=nokw,
    kw...
)
    filename = _maybe_add_suffix(filename, suffix)
    # This calls `create` in the /sources file for this `source`
    return create(filename, source, T, dims; lazy, missingval, kw...)
end
# Create from filename, type and extent with res or size keywords
function create(filename::Union{AbstractString,Nothing}, T::Union{Type,NamedTuple}, extent::Extents.Extent;
    res=nokw,
    size=nokw,
    crs=nothing,
    sampling=Points(),
    reverse_y=nokw,
    kw...
)
    dims = _extent2dims(extent; size, res, crs, sampling)
    dims = if reverse_y isa Bool && reverse_y && hasdim(ds, Y())
        DD.setdims(ds, reverse(DD.dims(ds, Y())))
    else
        dims
    end
    return create(filename, T, dims; kw...)
end
# Create in-memory Raster from type and dims
function create(filename::Nothing, ::Type{T}, dims::Tuple;
    missingval=nokw,
    fill=nokw,
    parent=nokw,
    verbose=true,
    f=identity,
    # Not used but here for consistency
    suffix=nokw,
    force=false,
    chunks=nokw,
    driver=nokw,
    options=nokw,
    kw...
) where T
    if verbose
        isnokw(chunks) || _warn_keyword_not_used("chunks", chunks)
        isnokw(driver) || _warn_keyword_not_used("driver", driver)
        isnokw(options) || _warn_keyword_not_used("options", options)
    end
    # Split inner and outer missingval if needed
    # For in-memory rasters we just ignore the inner value
    mv_inner, mv_outer = missingval isa Pair ? missingval : (missingval, missingval)
    # Get the element type from T and outer missingval
    eltype = isnokwornothing(mv_outer) ? T : promote_type(T, typeof(mv_outer))
    # Create the array
    data = if isnokwornothing(parent)
        Array{eltype}(undef, dims)
    else
        similar(parent, eltype, size(dims))
    end
    # Maybe fill the array
    isnokwornothing(fill) || fill!(data, fill)
    # Wrap as a Raster
    rast = Raster(data, dims; missingval=mv_outer, kw...)
    # Apply `f` before returning
    f(rast)
    return rast
end
# Create in-memory RasterStack from type and dims
function create(filename::Nothing, types::NamedTuple, dims::Tuple;
    suffix=keys(types),
    force=false,
    chunks=nokw,
    verbose=true,
    driver=nokw,
    options=nokw,
    parent=nokw,
    missingval=nokw,
    fill=nokw,
    layerdims=nokw,
    layermetadata=nokw,
    f=identity,
    kw...
)
    layerdims = isnokw(layerdims) ? map(_ -> basedims(dims), types) : layerdims
    layermetadata = layermetadata isa NamedTuple ? layermetadata : map(_ -> layermetadata, types)
    layerfill = fill isa NamedTuple ? fill : map(_ -> fill, types)
    layermissingvals = missingval isa NamedTuple ? missingval : map(_ -> missingval, types)
    layers = map(types, layermissingvals, layerfill, layerdims, layermetadata) do T, lmv, lfv, ld, lm
        create(nothing, T, DD.dims(dims, ld); 
            parent, missingval=lmv, fill=lfv, metadata=lm, driver, options,
        )
    end
    st = RasterStack(layers; kw...)
    f(st)
    return st
end
# Create on-disk Raster from filename, source, type and dims
function create(filename::AbstractString, source::Source, ::Type{T}, dims::Tuple;
    name=nokw,
    missingval=nokw,
    fill=nokw,
    metadata=nokw,
    chunks=nokw,
    scale=nokw,
    offset=nokw,
    dropband=!hasdim(dims, Band()),
    lazy=true,
    verbose=true,
    force=false,
    coerce=nokw,
    f=identity,
    kw...
) where T
    mv_inner, mv_outer = _missingval_pair(missingval)
    eltype = Missings.nonmissingtype(T)
    if isnokw(fill) || isnothing(fill)
        write = false # Leave fill undefined
        A = FillArrays.Zeros{eltype}(map(length, dims))
    else
        fill isa eltype || throw(ArgumentError("fill must be of type $eltype, got $fill"))
        write = true # Write fill to disk
        A = FillArrays.Fill{eltype}(fill, map(length, dims))
    end
    # Create layers of zero arrays
    rast = Raster(A, dims; name, missingval=mv_inner)
    # Write to disk
    Rasters.write(f, filename, source, rast;
        eltype, chunks, metadata, scale, offset, missingval=mv_inner, verbose, force, coerce, write, kw...
    ) do W
        # `write`` returns a variable, wrap it as a Raster and run f
        f(rebuild(rast; data=W))
    end
    # Don't pass in `missingval`, read it again from disk in case it changed
    return Raster(filename; source, lazy, metadata, dropband, coerce, missingval=mv_outer)
end
# Create on-disk RasterStack from filename, source, type and dims
function create(filename::AbstractString, source::Source, layertypes::NamedTuple, dims::Tuple;
    lazy=true,
    verbose=true,
    force=false,
    missingval=nokw,
    fill=nokw,
    metadata=nokw,
    layerdims=nokw,
    layermetadata=nokw,
    chunks=nokw,
    scale=nokw,
    offset=nokw,
    dropband=!hasdim(dims, Band),
    coerce=nokw,
    f=identity,
    kw...
)
    layerdims = if isnokwornothing(layerdims)
        # Use the same dims for all layers by default
        map(_ -> DD.basedims(dims), layertypes)
    else
        layerdims
    end
    mv_inner, mv_outer = _missingval_pair(missingval)
    # Only write value to disk variables if fill is defined
    write = !isnokwornothing(fill)
    # Make sure fill is per-layer
    fill_nt = fill isa NamedTuple ? fill : map(_ -> fill, layertypes)
    # Define no-allocation layers with FillArrays
    layers = map(layertypes, layerdims, fill_nt) do T, ld, fi
        lks = lookup(dims, ld)
        eltype = Missings.nonmissingtype(T)
        size = map(length, lks)
        if isnokwornothing(fi)
            A = FillArrays.Zeros{eltype}(size)
        else
            A = FillArrays.Fill{eltype}(fi, size)
        end
    end
    # Create layers of zero arrays
    st1 = RasterStack(layers, dims; layerdims, layermetadata, missingval=mv_inner)
    fn = Rasters.write(filename, st1;
        chunks, metadata, scale, offset, missingval=mv_inner, verbose, force, coerce, write, kw...
    ) do W
        # write returns a variable, wrap it as a RasterStack
        f(rebuild(st1; data=W))
    end
    st2 = RasterStack(fn; source, lazy, metadata, layerdims, dropband, coerce, missingval=mv_outer)
    return st2
end

_warn_keyword_not_used(label, obj) = @warn "`$label` of `$obj` found. But `chunks` are not used for in-memory rasters"
_missingval_pair(missingval::Pair) = missingval
_missingval_pair(missingval) = missingval => missingval