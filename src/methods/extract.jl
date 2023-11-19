"""
    extract(x, geoms; atol)

Extracts the value of `Raster` or `RasterStack` at given points, returning
an iterable of `NamedTuple` with properties for `:geometry` and raster or
stack layer values.

Note that if objects have more dimensions than the length of the point tuples,
sliced arrays or stacks will be returned instead of single values.

# Arguments

- `x`: a `Raster` or `RasterStack` to extract values from.
- `geoms`: GeoInterface.jl compatible geometries, or tables or iterables of geometries.

# Keywords

- `geometry`: include a `:geometry` column with the corresponding points for each value, `true` by default.
- `index`: include a column of the `CartesianIndex` for each value, `false` by default.
- `names`: `Tuple` of `Symbol` corresponding to layers of a `RasterStack`. All layers by default.
- `skipmissing`: skip missing points automatically.
- `atol`: a tolorerance for floating point lookup values for when the `LookupArray`
    contains `Points`. `atol` is ignored for `Intervals`.

# Example

Here we extact points matching the occurrence of the Mountain Pygmy Possum,
_Burramis parvus_. This could be used to fit a species distribution model.

```julia
using Rasters, RasterDataSources, ArchGDAL, GBIF2, CSV

# Get a stack of BioClim layers, and replace missing values with `missing`
st = RasterStack(WorldClim{BioClim}, (1, 3, 5, 7, 12)) |> replace_missing

# Download some occurrence data
obs = GBIF2.occurrence_search("Burramys parvus"; limit=5, year="2009")

# use `extract` to get values for all layers at each observation point.
# We `collect` to get a `Vector` from the lazy iterator.
collect(extract(st, pnts; skipmissing=true))

# output
5-element Vector{NamedTuple{(:geometry, :bio1, :bio3, :bio5, :bio7, :bio12)}}:
 (geometry = (0.21, 40.07), bio1 = 17.077084f0, bio3 = 41.20417f0, bio5 = 30.1f0, bio7 = 24.775f0, bio12 = 446.0f0)
 (geometry = (0.03, 39.97), bio1 = 17.076923f0, bio3 = 39.7983f0, bio5 = 29.638462f0, bio7 = 24.153847f0, bio12 = 441.0f0)
 (geometry = (0.03, 39.97), bio1 = 17.076923f0, bio3 = 39.7983f0, bio5 = 29.638462f0, bio7 = 24.153847f0, bio12 = 441.0f0)
 (geometry = (0.52, 40.37), bio1 = missing, bio3 = missing, bio5 = missing, bio7 = missing, bio12 = missing)
 (geometry = (0.32, 40.24), bio1 = 16.321388f0, bio3 = 41.659454f0, bio5 = 30.029825f0, bio7 = 25.544561f0, bio12 = 480.0f0)
```
"""
function extract end
function extract(x::RasterStackOrArray, data;
    dims=DD.dims(x, DEFAULT_POINT_ORDER), names=_names(x), kw...
)
    _extract(x, data; dims, names, kw...)
end
_extract(A::RasterStackOrArray, point::Missing; kw...) = missing
function _extract(A::RasterStackOrArray, geom; kw...)
    _extract(A, GI.geomtrait(geom), geom; kw...)
end
function _extract(A::RasterStackOrArray, ::Nothing, geoms; skipmissing=false, kw...)
    geom1 = first(Base.skipmissing(geoms))
    if GI.isgeometry(geom1) || GI.isfeature(geom1) || GI.isfeaturecollection(geom1)
        if skipmissing
            (_extract(A, g; kw...) for g in skipmissing(geoms))
        else
            (_extract(A, g; kw...) for g in geoms)
        end
    else
        throw(ArgumentError("`data` does not contain geomety objects"))
    end
end
function _extract(A::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature; kw...)
    _extract(A, GI.geometry(feature); kw...)
end
function _extract(A::RasterStackOrArray, ::GI.AbstractMultiPointTrait, geom; kw...)
    (_extract(A, p; kw...) for p in GI.getpoint(geom))
end
function _extract(A::RasterStackOrArray, ::GI.AbstractGeometryTrait, geom; 
    names, geometry=true, index=false, kw...
)
    B = boolmask(geom; to=dims(A, DEFAULT_POINT_ORDER), kw...)
    pts = DimPoints(B)
    dis = DimIndices(B)
    fs = let names=names
        ((b, i) -> _prop_nt(b, i, names),)
    end
    if geometry
        fs = (fs..., _geom_nt)
    end
    if index
        fs = (fs..., _index_nt)
    end
    return _create_nametuples(B, fs)
end
_geom_nt(x, pts) = NamedTuple{map(dim2key, dims(x))}(pts)
_prop_nt(st::AbstractRasterStack, I, names::NamedTuple{K}) where K = NamedTuple{K}(values(st[I]))
_prop_nt(A::AbstractRaster, I, names::NamedTuple{K}) where K = NamedTuple{K}((A[I],))
_index_nt(x, I) = (; index=CartesianIndex(I))

function _create_nametuple(B, fs::Tuple)
    (merge(map(f -> f(B, I), fs) for I in CartesianIndices(B) if B[I])
end

function _extract(x::RasterStackOrArray, ::GI.PointTrait, point; dims, names, atol=nothing)
    # Get the actual dimensions available in the object
    coords = map(DD.commondims(x, dims)) do d
        _dimcoord(d, point)
    end

    # Extract the values
    if any(map(ismissing, coords))
        # TODO test this branch somehow
        geometry = map(_ -> missing, coords)
        layer_vals = map(_ -> missing, names)
    else
        selectors = map(dims, coords) do d, c
            _at_or_contains(d, c, atol)
        end
        layer_vals = if DD.hasselection(x, selectors)
            x isa Raster ? (x[selectors...],) : x[selectors...]
        else
            map(_ -> missing, names)
        end
        geometry = point
    end
    properties = NamedTuple{keys(names)}(layer_vals)
    return (; geometry, properties...)
end

_names(A::AbstractRaster) = NamedTuple{(Symbol(name(A)),)}((Symbol(name(A)),))
_names(A::AbstractRasterStack) = NamedTuple{keys(A)}(keys(A))
