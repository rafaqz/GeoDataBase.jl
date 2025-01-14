# Object to hold extract keywords
struct Extractor{T,P,K,N<:NamedTuple{K},ID,G,I,S,F}
    names::N
    id::ID
    geometry::G
    index::I
    skipmissing::S
    flatten::F
end
Base.@constprop :aggressive @inline function Extractor(x, data;
    name::NTuple{<:Any,Symbol},
    skipmissing,
    flatten,
    geometry,
    id,
    index,
    kw...
)
    nt = NamedTuple{name}(name)
    id = _booltype(id)
    g = _booltype(geometry)
    i = _booltype(index)
    sm = _booltype(skipmissing)
    f = _booltype(flatten)
    Extractor(x, data, nt, id, g, i, sm, f)
end
function Extractor(
    x, data, nt::N, id::ID, g::G, i::I, sm::S, f::F
) where {N<:NamedTuple{K},ID,G,I,S,F} where K
    P = _proptype(x; skipmissing=sm, names=nt)
    T = _geom_rowtype(x, data; id, geometry=g, index=i, names=nt, skipmissing=sm)
    Extractor{T,P,K,N,ID,G,I,S,F}(nt, id, g, i, sm, f)
end

Base.eltype(::Extractor{T}) where T = T

# This mutable object is passed into closures as 
# fields are type-stable in clusores but variables are not
mutable struct LineRefs{T}
    i::Int
    prev::CartesianIndex{2}
    rows::Vector{T}
    function LineRefs{T}() where T
        i = 1
        prev = CartesianIndex((typemin(Int), typemin(Int)))
        new{T}(i, prev, Vector{T}())
    end
end

function _initialise!(lr::LineRefs{T}) where T
    lr.i = 1
    lr.prev = CartesianIndex((typemin(Int), typemin(Int)))
    lr.rows = Vector{T}()
end

_geom_rowtype(A, geom; kw...) =
    _geom_rowtype(A, GI.trait(geom), geom; kw...)
_geom_rowtype(A, ::Nothing, geoms; kw...) =
    _geom_rowtype(A, first(skipmissing(geoms)); kw...)
_geom_rowtype(A, ::GI.AbstractGeometryTrait, geom; kw...) =
    _geom_rowtype(A, first(GI.getpoint(geom)); kw...)
_geom_rowtype(A, ::Nothing, geoms::Missing; kw...) =
    _rowtype(A, missing; kw...)
function _geom_rowtype(A, ::GI.AbstractPointTrait, p; kw...)
    tuplepoint = if GI.is3d(p)
        (GI.x(p), GI.y(p), GI.z(p))
    else
        (GI.x(p), GI.y(p))
    end
    return _rowtype(A, tuplepoint; kw...)
end

# skipinvalid: can G and I be missing. skipmissing: can nametypes be missing
function _proptype(x;
    skipmissing, names::NamedTuple{K}, kw...
) where K
    NamedTuple{K,Tuple{_nametypes(x, names, skipmissing)...}}
end

"""
    extract(x, geometries; kw...)

Extracts the value of `Raster` or `RasterStack` for the passed in geometries,
returning an `Vector{NamedTuple}` with properties for `:geometry` and `Raster`
or `RasterStack` layer values.

For lines, linestrings and linear rings points are extracted for each pixel that
the line touches.

For polygons, all cells witih centers covered by the polygon are returned.

Note that if objects have more dimensions than the length of the point tuples,
sliced arrays or stacks will be returned instead of single values.

# Arguments

- `x`: a `Raster` or `RasterStack` to extract values from.
$DATA_ARGUMENT

# Keywords

- `geometry`: include a `:geometry` field in rows, which will be a
    tuple point. Either the original point for points or the pixel
    center point for line and polygon extract. `true` by default.
- `index`: include `:index` field of extracted points in rows, `false` by default.
- `name`: a `Symbol` or `Tuple` of `Symbol` corresponding to layer/s
    of a `RasterStack` to extract. All layers are extracted by default.
- `skipmissing`: skip missing points automatically.
- `flatten`: flatten extracted points from multiple 
    geometries into a single vector. `true` by default. 
    Unmixed point geometries are always flattened.
    Flattening is slow and single threaded, `flatten=false` may be a 
    large performance improvement in combination with `threaded=true`.
- `atol`: a tolerance for floating point lookup values for when the `Lookup`
    contains `Points`. `atol` is ignored for `Intervals`.
$BOUNDARY_KEYWORD
$GEOMETRYCOLUMN_KEYWORD

# Example

Here we extract points matching the occurrence of the Mountain Pygmy Possum,
_Burramis parvus_. This could be used to fit a species distribution model.

```julia
using Rasters, RasterDataSources, ArchGDAL, GBIF2, CSV

# Get a stack of BioClim layers, and replace missing values with `missing`
st = RasterStack(WorldClim{BioClim}, (1, 3, 5, 7, 12)) |> replace_missing

# Download some occurrence data
obs = GBIF2.occurrence_search("Burramys parvus"; limit=5, year="2009")

# use `extract` to get values for all layers at each observation point.
# We `collect` to get a `Vector` from the lazy iterator.
extract(st, obs; skipmissing=true)

# output
5-element Vector{NamedTuple{(:geometry, :bio1, :bio3, :bio5, :bio7, :bio12)}}:
 (geometry = (0.21, 40.07), bio1 = 17.077084f0, bio3 = 41.20417f0, bio5 = 30.1f0, bio7 = 24.775f0, bio12 = 446.0f0)
 (geometry = (0.03, 39.97), bio1 = 17.076923f0, bio3 = 39.7983f0, bio5 = 29.638462f0, bio7 = 24.153847f0, bio12 = 441.0f0)
 (geometry = (0.03, 39.97), bio1 = 17.076923f0, bio3 = 39.7983f0, bio5 = 29.638462f0, bio7 = 24.153847f0, bio12 = 441.0f0)
 (geometry = (0.52, 40.37), bio1 = missing, bio3 = missing, bio5 = missing, bio7 = missing, bio12 = missing)
 (geometry = (0.32, 40.24), bio1 = 16.321388f0, bio3 = 41.659454f0, bio5 = 30.029825f0, bio7 = 25.544561f0, bio12 = 480.0f0)
```

Note: passing in arrays, geometry collections or feature collections
containing a mix of points and other geometries has undefined results.
"""
function extract end
Base.@constprop :aggressive @inline function extract(x::RasterStackOrArray, data;
    geometrycolumn=nothing,
    names::NTuple{<:Any,Symbol}=_names(x),
    name::NTuple{<:Any,Symbol}=names,
    skipmissing=false,
    flatten=true,
    id=false,
    geometry=true,
    index=false,
    kw...
)
    n = DD._astuple(name)
    g, g1 = if GI.isgeometry(data)
        data, data
    elseif GI.isfeature(data)
        g = GI.geometry(data)
        g, g
    else
        gs = _get_geometries(data, geometrycolumn)
        gs, (isempty(gs) || all(ismissing, gs)) ? missing : first(Base.skipmissing(gs))
    end

    xp = _prepare_for_burning(x)
    e = Extractor(xp, g1; name, skipmissing, flatten, id, geometry, index)
    id_init = 1
    return _extract(xp, e, id_init, g; kw...)
end

# TODO use a GeometryOpsCore method like `applyreduce` here?
function _extract(A::RasterStackOrArray, e::Extractor{T}, id::Int, geom::Missing; kw...) where T
    return if istrue(e.skipmissing)
        T[]
    else
        T[_maybe_add_fields(T, map(_ -> missing, e.names), id, missing, missing)]
    end
end
function _extract(A::RasterStackOrArray, e::Extractor, id::Int, geom; kw...)
    _extract(A, e, id, GI.geomtrait(geom), geom; kw...)
end
function _extract(A::RasterStackOrArray, e::Extractor{T}, id::Int, ::Nothing, geoms::AbstractArray;
    threaded=false, progress=true, kw...
) where T
    # Handle empty / all missing cases
    isempty(geoms) && return T[]
    geom1 = all(ismissing, geoms) ? missing : first(Base.skipmissing(geoms))
    trait1 = GI.trait(geom1)
    # We split out points from other geoms for performance
    if trait1 isa GI.PointTrait
        allpoints = true
        rows = Vector{T}(undef, length(geoms))
        if istrue(e.skipmissing)
            if threaded
                nonmissing = Vector{Bool}(undef, length(geoms))
                _run(1:length(geoms), threaded, progress, "Extracting points...") do i
                    g = geoms[i]
                    GI.geomtrait(g) isa GI.PointTrait || return nothing
                    loc_id = id + i - 1
                    nonmissing[i] = _extract_point!(rows, A, e, loc_id, g, i; kw...)::Bool
                    return nothing
                end
                # This could use less memory if we reuse `rows` and shorten it
                rows = rows[nonmissing]
            else
                j_ref = Ref(1)
                # For non-threaded be more memory efficient
                _run(1:length(geoms), false, progress, "Extracting points...") do i
                    g = geoms[i]
                    loc_id = id + i - 1
                    if GI.geomtrait(g) isa GI.PointTrait
                        j_ref[] += _extract_point!(rows, A, e, loc_id, g, j_ref[]; kw...)::Bool
                    end
                    return nothing
                end
                deleteat!(rows, j_ref[]:length(rows))
            end
        else
            _run(1:length(geoms), threaded, progress, "Extracting points...") do i
                g = geoms[i]
                loc_id = id + i - 1
                _extract_point!(rows, A, e, loc_id, g, i; kw...)
                return nothing
            end
        end
        # If we found a non-point geometry,
        # ignore these rows and start again generically
        allpoints && return rows
    end
    
    row_vecs = Vector{Vector{T}}(undef, length(geoms))
    if threaded
        thread_line_refs = [LineRefs{T}() for _ in 1:Threads.nthreads()]
        _run(1:length(geoms), threaded, progress, "Extracting geometries...") do i
            line_refs = thread_line_refs[Threads.threadid()]
            loc_id = id + i - 1
            row_vecs[i] = _extract(A, e, loc_id, geoms[i]; line_refs, kw...)
        end
    else
        line_refs = LineRefs{T}()
        _run(1:length(geoms), threaded, progress, "Extracting geometries...") do i
            loc_id = id + i - 1
            row_vecs[i] = _extract(A, e, loc_id, geoms[i]; line_refs, kw...)
        end
    end

    # TODO this is a bit slow and only on one thread
    if istrue(e.flatten)
        n = sum(map(length, row_vecs))
        out_rows = Vector{T}(undef, n)
        k::Int = 1
        for row_vec in row_vecs
            for row in row_vec
                out_rows[k] = row
                k += 1
            end
        end
        return out_rows
    else
        return row_vecs
    end
end
function _extract(
    A::RasterStackOrArray, e::Extractor{T}, id::Int, ::GI.AbstractMultiPointTrait, geom; kw...
)::Vector{T} where T
    n = GI.npoint(geom)
    rows = _init_rows(e, n)
    i = 1
    for p in GI.getpoint(geom)
        i += _extract_point!(rows, A, e, id, p, i; kw...)::Bool
    end
    if istrue(skipmissing)
        # Remove excees rows where missing
        deleteat!(rows, i:length(rows))
    end
    return rows
end
function _extract(A::RasterStackOrArray, e::Extractor{T}, id::Int, ::GI.PointTrait, p; kw...) where T
    rows = _init_rows(e, 1)
    _extract_point!(rows, A, e, id, p, 1; kw...)
end
@noinline function _extract(
    A::RasterStackOrArray, e::Extractor{T}, id::Int, ::GI.AbstractLineStringTrait, geom;
    line_refs=LineRefs{T}(),
    kw...
)::Vector{T} where T
    _initialise!(line_refs)
    # Subset/offst is not really needed for line buring
    # But without it we can get different fp errors
    # to polygons and eng up with lines in different
    # places when they are right on the line.
    dims, offset = _template(A, geom)
    dp = DimPoints(A)
    function _length_callback(n)
        resize!(line_refs.rows, n)
    end

    _burn_lines!(_length_callback, dims, geom) do D
        I = CartesianIndex(map(val, D))
        # Avoid duplicates from adjacent line segments
        line_refs.prev == I && return nothing
        line_refs.prev = I
        # Make sure we only hit this pixel once
        # D is always inbounds
        line_refs.i += _maybe_set_row!(line_refs.rows, e.skipmissing, e, id, A, dp, offset, I, line_refs.i)
        return nothing
    end
    deleteat!(line_refs.rows, line_refs.i:length(line_refs.rows))
    return line_refs.rows
end
@noinline function _extract(
    A::RasterStackOrArray, e::Extractor{T}, id::Int, ::GI.AbstractGeometryTrait, geom; kw...
)::Vector{T} where T
    # Make a raster mask of the geometry
    dims, offset = _template(A, geom)
    B = boolmask(geom; to=dims, burncheck_metadata=NoMetadata(), kw...)
    n = count(B)
    dp = DimPoints(A)
    i = 1
    rows = _init_rows(e, n)
    # Add a row for each pixel that is `true` in the mask
    for I in CartesianIndices(B)
        B[I] || continue
        i += _maybe_set_row!(rows, e.skipmissing, e, id, A, dp, offset, I, i)
    end
    # Cleanup
    deleteat!(rows, i:length(rows))
    return rows
end

function _template(x, geom)
    ods = otherdims(x, DEFAULT_POINT_ORDER)
    # Build a dummy raster in case this is a stack
    # views are just easier on an array than dims
    t1 = Raster(FillArrays.Zeros(size(x)), dims(x))
    t2 = if length(ods) > 0
        view(t1, map(d -> rebuild(d, firstindex(d)), ods))
    else
        t1
    end
    I = dims2indices(dims(t2), Touches(GI.extent(geom)))
    t3 = view(t2, I...)
    offset = CartesianIndex(map(i -> first(i) - 1, I))
    return dims(t3), offset
end

Base.@assume_effects :foldable function _maybe_set_row!(
    rows, skipmissing::_True, e, id, A, dp, offset, I, i;
)
    props = _prop_nt(e, A, I, skipmissing)
    return if ismissing(props)
        0
    else
        _maybe_set_row!(rows, _False(), e, id, A, dp, offset, I, i; props)
    end
end
Base.@assume_effects :foldable function _maybe_set_row!(
    rows::Vector{T}, skipmissing::_False, e, id, A, dp, offset, I, i;
    props=_prop_nt(e, A, I, skipmissing)
) where T
    Ioff = I + offset
    geom = dp[Ioff]
    rows[i] = _maybe_add_fields(T, props, id, geom, Tuple(Ioff))
    return 1
end

_init_rows(e::Extractor{T}, n=0) where T = Vector{T}(undef, n)

Base.@assume_effects :foldable _ismissingval(A::Union{Raster,RasterStack}, props)::Bool =
    _ismissingval(missingval(A), props)
Base.@assume_effects :foldable _ismissingval(A::Union{Raster,RasterStack}, props::NamedTuple)::Bool =
    _ismissingval(missingval(A), props)
Base.@assume_effects :foldable _ismissingval(mvs::NamedTuple{K}, props::NamedTuple{K}) where K =
    any(DD.unrolled_map((p, mv) -> ismissing(p) || p === mv, props, mvs))::Bool
Base.@assume_effects :foldable _ismissingval(mvs::NamedTuple, props::NamedTuple{K}) where K =
    _ismissingval(mvs[K], props)::Bool
Base.@assume_effects :foldable _ismissingval(mv, props::NamedTuple)::Bool =
    any(DD.unrolled_map(p -> ismissing(p) || p === mv, props))
Base.@assume_effects :foldable _ismissingval(mv, prop)::Bool = (mv === prop)

# We always return NamedTuple
Base.@assume_effects :foldable function _prop_nt(
    ::Extractor{<:Any,P,K}, st::AbstractRasterStack, I, sm::_False
)::P where {P,K}
    P(st[K][I])
end
Base.@assume_effects :foldable function _prop_nt(
    ::Extractor{<:Any,P,K}, st::AbstractRasterStack{K}, I, sm::_False
)::P where {P,K}
    P(st[I])
end
Base.@assume_effects :foldable function _prop_nt(
    ::Extractor{<:Any,P,K}, A::AbstractRaster, I, sm::_False
)::P where {P,K}
    P(NamedTuple{K}((A[I],)))
end
Base.@assume_effects :foldable function _prop_nt(
    ::Extractor{<:Any,P,K}, st::AbstractRasterStack, I, sm::_True
)::Union{P,Missing} where {P,K}
    x = st[K][I]
    _ismissingval(A, x) ? missing : x::P
end
Base.@assume_effects :foldable function _prop_nt(
    ::Extractor{<:Any,P,K}, st::AbstractRasterStack{K}, I, sm::_True
)::Union{P,Missing} where {P,K}
    x = st[I]
    _ismissingval(st, x) ? missing : x::P
end
Base.@assume_effects :foldable function _prop_nt(
    ::Extractor{<:Any,P,K}, A::AbstractRaster, I, sm::_True
)::Union{P,Missing} where {P,K}
    x = A[I]
    _ismissingval(A, x) ? missing : NamedTuple{K}((x,))::P
end

# Extract a single point
# Missing point to remove
@inline function _extract_point!(rows::Vector{T}, x::RasterStackOrArray, e::Extractor, id, point, i; kw...) where T
    _extract_point!(rows, x, e, id, e.skipmissing, point, i; kw...)
end
@inline function _extract_point!(
    rows::Vector{T}, x::RasterStackOrArray, e::Extractor{T,<:Any,K}, id, skipmissing::_True, point::Missing, i;
    kw...
) where {T,K}
    return false
end
# Normal point with missing / out of bounds data removed
@inline function _extract_point!(
    rows::Vector{T}, x::RasterStackOrArray, e::Extractor{T,<:Any,K}, id, skipmissing::_True, point, i;
    dims=DD.dims(x, DEFAULT_POINT_ORDER),
    atol=nothing,
    kw...
) where {T,K}
    # Get the actual dimensions available in the object
    coords = map(DD.dims(x, dims)) do d
        _dimcoord(d, point)
    end
    # If any are coordinates missing, also return missing for everything
    if !any(map(ismissing, coords))
        selector_dims = map(dims, coords) do d, c
            _at_or_contains(d, c, atol)
        end
        selectors = map(val, DD.dims(selector_dims, DD.dims(x)))
        # Check the selector is in bounds / actually selectable
        I = DD.selectindices(DD.dims(x, dims), selectors; err=_False())::Union{Nothing,Tuple{Int,Int}}
        if !isnothing(I)
            D = map(rebuild, selector_dims, I)
            if x isa Raster
                prop = x[D]
                _ismissingval(missingval(x), prop) && return false
                props = NamedTuple{K,Tuple{eltype(x)}}((prop,))
            else
                props = x[D][K]
                _ismissingval(missingval(x), props) && return false
            end
            rows[i] = _maybe_add_fields(T, props, id, coords, I)
            return true
        end
    end
    return false
end
# Missing point to keep
@inline function _extract_point!(
    rows::Vector{T}, x::RasterStackOrArray, e::Extractor{T,P,K}, id, skipmissing::_False, point::Missing, i; kw...
) where {T,P,K}
    props = map(_ -> missing, e.names)
    geom = missing
    I = missing
    rows[i] = _maybe_add_fields(T, props, id, geom, I)
    return true
end
# Normal point with missing / out of bounds data kept with `missing` fields
@inline function _extract_point!(
    rows::Vector{T}, x::RasterStackOrArray, e::Extractor{T,P,K}, id, skipmissing::_False, point, i;
    dims=DD.dims(x, DEFAULT_POINT_ORDER),
    atol=nothing,
    kw...
) where {T,P,K}
    # Get the actual dimensions available in the object
    coords = map(dims) do d
        _dimcoord(d, point)
    end
    # If any are coordinates missing, also return missing for everything
    rows[i] = if any(map(ismissing, coords))
        _maybe_add_fields(T, map(_ -> missing, names), id, missing, missing)
    else
        selector_dims = map(dims, coords) do d, c
            _at_or_contains(d, c, atol)
        end
        selectors = map(val, DD.dims(selector_dims, DD.dims(x)))
        # Check the selector is in bounds / actually selectable
        I = DD.selectindices(DD.dims(x, dims), selectors; err=_False())::Union{Nothing,Tuple{Int,Int}}
        if isnothing(I)
            _maybe_add_fields(T, map(_ -> missing, e.names), id, coords, missing)
        else
            D = map(rebuild, selector_dims, I)
            props = if x isa Raster
                P(NamedTuple{K,Tuple{eltype(x)}}((x[D],)))
            else
                P(NamedTuple(x[D])[K])
            end
            _maybe_add_fields(T, props, id, coords, I)
        end
    end
    return 1
end

# Maybe add optional fields
# It is critically important for performance that this is type stable
Base.@assume_effects :total function _maybe_add_fields(
    ::Type{T}, props::NamedTuple, id, point::Union{Tuple,Missing}, I
)::T where {T<:NamedTuple{K}} where K
    row = :index in K ? merge((; index=I), props) : props 
    row = :geometry in K ? merge((; geometry=point), row) : row
    row = :id in K ? merge((; id), row) : row
    return T(row)
end
