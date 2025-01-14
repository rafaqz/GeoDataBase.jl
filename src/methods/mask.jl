
const INVERT_KEYWORD = """
- `invert`: invert the mask, so that areas no missing in `with` are
    masked, and areas missing in `with` are masked.
"""

const COLLAPSE_KEYWORD = """
- `collapse`: if `true`, collapse all geometry masks into a single mask. Otherwise
    return a Raster with an additional `geometry` dimension, so that each slice
    along this axis is the mask of the `geometry` opbject of each row of the
    table, feature in the feature collection, or just each geometry in the iterable.
"""

"""
    mask(A:AbstractRaster; with, missingval=missingval(A))
    mask(x; with)

Return a new array with values of `A` masked by the missing values of `with`,
or by the shape of `with`, if `with` is a geometric object.

# Arguments

- `x`: a `Raster` or `RasterStack`

# Keywords

- `with`: an `AbstractRaster`, or any GeoInterface.jl compatible objects
    or table. The coordinate reference system of the point must match `crs(A)`.
$INVERT_KEYWORD
- `missingval`: the missing value to use in the returned file.
$FILENAME_KEYWORD
$SUFFIX_KEYWORD

# Geometry keywords

These can be used when `with` is a GeoInterface.jl compatible object:

$SHAPE_KEYWORD
$BOUNDARY_KEYWORD
$GEOMETRYCOLUMN_KEYWORD

# Example

Mask an unmasked AWAP layer with a masked WorldClim layer,
by first resampling the mask.

```julia
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates

# Load and plot the file
awap = read(Raster(AWAP, :tmax; date=DateTime(2001, 1, 1)))
a = plot(awap; clims=(10, 45))

# Create a mask my resampling a worldclim file
wc = Raster(WorldClim{Climate}, :prec; month=1)
wc_mask = resample(wc; to=awap)

# Mask
awap_masked = mask(awap; with=wc_mask)
b = plot(awap_masked; clims=(10, 45))

savefig(a, "build/mask_example_before.png");
savefig(b, "build/mask_example_after.png"); nothing
# output

```

### Before `mask`:


### After `mask`:


$EXPERIMENTAL
"""
mask(x; with, kw...) = _mask(x, with; kw...)

# Geometry mask
function _mask(s::AbstractRasterSeries, with; kw...)
    B = boolmask(with; to=dims(first(s), DEFAULT_POINT_ORDER), kw...)
    return _mask(s, B)
end
function _mask(x::RasterStackOrArray, with; kw...)
    # Geometries can only have `X`/`Y`/`Z` dims so limit them here
    B = boolmask(with; to=dims(x, DEFAULT_POINT_ORDER), kw...)
    return _mask(x, B)
end
# Array mask
function _mask(A::AbstractRaster, with::AbstractRaster;
    filename=nothing, suffix=nothing, missingval=_missingval_or_missing(A), kw...
)
    missingval = ismissing(missingval) ? missing : convert(eltype(A), missingval)
    return create(filename, A; suffix, missingval) do C
        # The values array will be be written to A1 in `mask!`
        mask!(C; with, missingval, values=A, kw...)
    end
end
function _mask(xs::AbstractRasterStack, with::AbstractRaster; suffix=keys(xs), kw...)
    mapargs((x, s) -> mask(x; with, suffix=s, kw...), xs, suffix)
end
function _mask(xs::AbstractRasterSeries, with::AbstractRaster; kw...)
    map(x -> mask(x; with, kw...), xs)
end

"""
    mask!(x; with, missingval=missingval(A))

Mask `A` by the missing values of `with`, or by all values outside `with` if it is a polygon.

If `with` is a polygon, creates a new array where points falling outside the polygon
have been replaced by `missingval(A)`.

Return a new array with values of `A` masked by the missing values of `with`,
or by a polygon.

# Arguments

- `x`: a `Raster` or `RasterStack`.

# Keywords

- `with`: another `AbstractRaster`, a `AbstractVector` of `Tuple` points,
    or any GeoInterface.jl `AbstractGeometry`. The coordinate reference system
    of the point must match `crs(A)`.
$INVERT_KEYWORD
- `missingval`: the missing value to write to A in masked areas,
    by default `missingval(A)`.

# Example

Mask an unmasked AWAP layer with a masked WorldClim layer,
by first resampling the mask to match the size and projection.

```julia
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates

# Load and plot the file
awap = read(RasterStack(AWAP, (:tmin, :tmax); date=DateTime(2001, 1, 1)))
a = plot(awap; clims=(10, 45), c=:imola)

# Create a mask my resampling a worldclim file
wc = Raster(WorldClim{Climate}, :prec; month=1)
wc_mask = resample(wc; to=awap)

# Mask
mask!(awap; with=wc_mask)
b = plot(awap; clims=(10, 45))

savefig(a, "build/mask_bang_example_before.png");
savefig(b, "build/mask_bang_example_after.png"); nothing

# output

```

### Before `mask!`:


### After `mask!`:


$EXPERIMENTAL
"""
function mask! end
function mask!(xs::AbstractRasterSeries; kw...)
    foreach(x -> mask!(x; kw...), xs)
    return xs
end
mask!(xs::RasterStackOrArray; with, kw...) = _mask!(xs, with; kw...)

# Geometry mask
function _mask!(x::RasterStackOrArray, geom; kw...)
    B = boolmask(geom; to=dims(x), kw...)
    _mask!(x, B; missingval=false)
    return x
end
# Array mask
function _mask!(st::RasterStack, with::AbstractRaster; kw...)
    maplayers(A -> mask!(A; with, kw...), st)
    return st
end

function _mask!(A::AbstractRaster, with::AbstractRaster; 
    missingval=missingval(A), 
    values=A,
    invert=false,
)
    missingval isa Nothing && _nomissingerror()
    missingval = convert(eltype(A), missingval)

    if invert
        broadcast_dims!(A, values, with) do x, w
            if (ismissing(w) || isequal(w, Rasters.missingval(with))) && 
                !(ismissing(x) || isequal(x, Rasters.missingval(values)))
                convert(eltype(A), x)
            else
                missingval
            end
        end
    else
        broadcast_dims!(A, values, with) do x, w
            if (ismissing(w) || isequal(w, Rasters.missingval(with))) || 
                (ismissing(x) || isequal(x, Rasters.missingval(values)))
                missingval
            else
                convert(eltype(A), x)
            end
        end
    end
    return rebuild(A, missingval=missingval)
end

_nomissingerror() = throw(ArgumentError("Array has no `missingval`. Pass a `missingval` keyword compatible with the type, or use `rebuild(A; missingval=somemissingval)` to set it."))

"""
    boolmask(obj::Raster; [missingval])
    boolmask(obj; [to, res, size])
    boolmask(obj::RasterStack; alllayers=true, kw...)

Create a mask array of `Bool` values, from another `Raster`.
`AbstractRasterStack` or `AbstractRasterSeries` are also accepted. 

The array returned from calling `boolmask` on a `AbstractRaster` is a
[`Raster`](@ref) with the same dimensions as the original array and a
`missingval` of `false`.

# Arguments

- $OBJ_ARGUMENT

# `Raster` / `RasterStack` Keywords

$INVERT_KEYWORD
- `missingval`: The missing value of the source array, with default `missingval(raster)`.

# Keywords

- `alllayers`: if `true` a mask is taken for all layers, otherwise only the first layer is used. Defaults to `true`

$GEOM_KEYWORDS
$GEOMETRYCOLUMN_KEYWORD
$THREADED_KEYWORD
$PROGRESS_KEYWORD

For tabular data, feature collections and other iterables

$COLLAPSE_KEYWORD

# Example

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates
wc = Raster(WorldClim{Climate}, :prec; month=1)
boolmask(wc) |> plot

savefig("build/boolmask_example.png"); nothing

# output
```

![boolmask](boolmask_example.png)

$EXPERIMENTAL
"""
function boolmask end

function boolmask(stack::AbstractRasterStack; 
    alllayers=true, 
    invert=false,
    to=dims(stack), 
    kw...
) 
    if alllayers
        _mask_multilayer(stack, to; 
            invert, kw..., _dest_presentval=!invert, _dest_missingval=invert
        )
    else
        boolmask(layers(stack, 1); invert, kw...)
    end
end

function boolmask(series::AbstractRasterSeries; 
    alllayers=true, 
    to=first(series), 
    invert=false,
    kw...
)
    if alllayers
        _mask_multilayer(series, to; 
            invert, kw..., _dest_presentval=!invert, _dest_missingval=invert
        )
    else
        boolmask(first(series); invert, kw...)
    end
end

function boolmask(source::AbstractRaster; 
    invert::Bool=false, 
    kw...
)
    dest = _init_bools(source, BitArray, nothing; kw..., missingval=invert)
    boolmask!(dest, source; invert, kw...)
    return rebuild(dest; missingval=false)
end
# this method is used where x is a geometry
function boolmask(x; 
    to=nothing, 
    invert::Bool=false,
    geometrycolumn=nothing,
    kw...
)
    if to isa Union{AbstractDimArray,AbstractDimStack,DimTuple}
        to = dims(to, DEFAULT_POINT_ORDER)
    end
    data = isnothing(GI.geomtrait(x)) ? _get_geometries(x, geometrycolumn) : x
    dest = _init_bools(to, BitArray, data; kw..., missingval=invert)
    boolmask!(dest, data; invert, kw...)
    return rebuild(dest; missingval=false)
end

function boolmask!(dest::AbstractRaster, src::AbstractRaster;
    missingval=_missingval_or_missing(src),
    invert=false,
)
    if invert
        broadcast_dims!(x -> isequal(x, missingval), dest, src)
    else
        broadcast_dims!(x -> !isequal(x, missingval), dest, src)
    end
end
function boolmask!(dest::AbstractRaster, data;
    invert=false,
    lock=nothing, 
    progress=true, 
    threaded=false, 
    burncheck_metadata=Metadata(),
    allocs=_burning_allocs(dest; threaded, burncheck_metadata), 
    geometrycolumn=nothing,
    collapse=nokw,
    kw...
)
    if collapse === false && hasdim(dest, :geometry)
        geoms = _get_geometries(data, geometrycolumn)
        range = eachindex(geoms)
        _run(range, threaded, progress, "Burning each geometry to a BitArray slice...") do i
            geom = geoms[i]
            ismissing(geom) && return nothing
            slice = view(dest, Dim{:geometry}(i))
            # We don't need locks - these are independent slices
            burn_geometry!(slice, geom; kw..., fill=!invert, allocs=_get_alloc(allocs))
            return nothing
        end
    elseif isnokw(collapse) || collapse === true
        burn_geometry!(dest, data; kw..., allocs, lock, progress, threaded, geometrycolumn, fill=!invert)
    else
        throw(ArgumentError("`collapse` must be `false` or not passed if there is no `:geometry` dimension in `dest`"))
    end
    return dest
end

"""
    missingmask(obj::Raster; kw...)
    missingmask(obj; [to, res, size])
    missingmask(obj::RasterStack; alllayers=true, kw...)

Create a mask array of `missing` and `true` values, from another `Raster`.
`AbstractRasterStack` or `AbstractRasterSeries` are also accepted-

For [`AbstractRaster`](@ref) the default `missingval` is `missingval(A)`,
but others can be chosen manually.

The array returned from calling `missingmask` on a `AbstractRaster` is a
[`Raster`](@ref) with the same size and fields as the original array.

# Arguments

- `obj`: $OBJ_ARGUMENT

# Keywords

- `alllayers`: if `true` a mask is taken for all layers, otherwise only the first layer is used. Defaults to `true`
$INVERT_KEYWORD
$GEOM_KEYWORDS
$GEOMETRYCOLUMN_KEYWORD

# Example

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates
wc = Raster(WorldClim{Climate}, :prec; month=1)
missingmask(wc) |> plot

savefig("build/missingmask_example.png"); nothing

# output
```

![missingmask](missingmask_example.png)

$EXPERIMENTAL
"""
function missingmask(x; kw...)
    B = boolmask(x; kw...)
    M = Array{Union{Missing,Bool}}(undef, size(B))
    M .= _false_to_missing.(B)
    return rebuild(B; data=M, missingval=missing)
end
function missingmask(source::AbstractRaster; kw...)
    dest = _init_bools(source, Array{Union{Missing,Bool}}, nothing; kw..., missingval=missing)
    return rebuild(missingmask!(dest, source; kw...); missingval=missing)
end

function missingmask!(dest::AbstractRaster, src::AbstractRaster;
    invert=false,
    missingval=_missingval_or_missing(src)
)
    if invert
        broadcast_dims!(dest, src) do x
            isequal(x, missingval) ? true : missing
        end
    else
        broadcast_dims!(dest, src) do x
            isequal(x, missingval) ? missing : true
        end
    end
end
function missingmask!(dest::AbstractRaster, geom; kw...)
    # boolmask! handles `invert` keyword here
    B = boolmask!(dest, geom; kw...) 
    dest .= _false_to_missing.(B)
    return dest
end

_false_to_missing(b::Bool) = (b ? true : missing)::Union{Missing,Bool}

_mask_multilayer(st::AbstractRasterStack, to; kw...) =
    _mask_multilayer(layers(st), to; kw...)
function _mask_multilayer(layers::Union{<:NamedTuple,<:AbstractRasterSeries}, to; 
    _dest_presentval,
    _dest_missingval, 
    missingval=nothing, 
    kw...
)
    T = Union{typeof(_dest_missingval),Bool}
    dest = _init_bools(to, Array{T}, layers; kw..., missingval=_dest_missingval)
    dest .= _dest_presentval

    missingval = if isnothing(missingval)
        map(_missingval_or_missing, layers)
    elseif missingval isa NamedTuple
        missingval
    else
        map(_ -> missingval, layers)
    end

    map(layers, missingval) do layer, mv
        broadcast_dims!(dest, dest, layer) do d, x
            isequal(d, _dest_missingval) || isequal(x, mv) ? _dest_missingval : _dest_presentval
        end
    end
   return dest
end
