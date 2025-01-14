"""
    AbstractProjected <: AbstractSampled

Abstract supertype for projected index lookups.
"""
abstract type AbstractProjected{T,O,Sp,Sa} <: AbstractSampled{T,O,Sp,Sa} end

struct AutoDim end

# For now we just remove CRS on GPU - it often contains strings
function Adapt.adapt_structure(to, l::AbstractProjected)
    sampled = Sampled(parent(l), order(l), span(l), sampling(l), metadata(l))
    return Adapt.adapt_structure(to, sampled)
end

GeoInterface.crs(lookup::Lookup) = nothing
mappedcrs(lookup::Lookup) = nothing

# When the lookup is formatted with an array we match the `dim` field with the
# wrapper dimension. We will need this later for e.g. projecting with GDAL,
# where we need to know which lookup is X and which is Y
function Dimensions.format(m::AbstractProjected, D::Type, index, axis::AbstractRange)
    i = Dimensions._format(index, axis)
    o = Dimensions._format(order(m), D, index)
    sp = Dimensions._format(span(m), D, index)
    sa = Dimensions._format(sampling(m), sp, D, index)
    dim = Dimensions.basetypeof(D)()
    x = rebuild(m; data=i, order=o, span=sp, sampling=sa, dim=dim)
    return x
end

"""
    Projected <: AbstractProjected

    Projected(order, span, sampling, crs, mappedcrs)
    Projected(; order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(), crs, mappedcrs=nothing)

An [`AbstractSampled`]($DDabssampleddocs) `Lookup` with projections attached.

Fields and behaviours are identical to [`Sampled`]($DDsampleddocs)
with the addition of `crs` and `mappedcrs` fields.

If both `crs` and `mappedcrs` fields contain CRS data (in a `GeoFormat` wrapper
from GeoFormatTypes.jl) the selector inputs and plot axes will be converted
from and to the specified `mappedcrs` projection automatically. A common use case
would be to pass `mappedcrs=EPSG(4326)` to the constructor when loading eg. a GDALarray:

```julia
GDALarray(filename; mappedcrs=EPSG(4326))
```

The underlying `crs` will be detected by GDAL.

If `mappedcrs` is not supplied (ie. `mappedcrs=nothing`), the base index will be
shown on plots, and selectors will need to use whatever format it is in.
"""
struct Projected{T,A<:AbstractVector{T},O<:Order,Sp<:Span,Sa<:Sampling,MD,PC<:Union{GeoFormat,Nothing},MC<:Union{GeoFormat,Nothing},D} <: AbstractProjected{T,O,Sp,Sa}
    data::A
    order::O
    span::Sp
    sampling::Sa
    metadata::MD
    crs::PC
    mappedcrs::MC
    dim::D
end
function Projected(data=AutoValues();
    order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(),
    metadata=NoMetadata(),
    crs::Union{GeoFormat,Nothing},
    mappedcrs::Union{GeoFormat,Nothing}=nothing,
    dim=AutoDim()
)
    Projected(data, order, span, sampling, metadata, crs, mappedcrs, dim)
end
function Projected(l::Sampled;
    order=order(l), span=span(l), sampling=sampling(l),
    metadata=metadata(l),
    crs::Union{GeoFormat,Nothing,NoKW},
    mappedcrs::Union{GeoFormat,NoKW,Nothing}=nokw,
    dim=AutoDim()
)
    crs = isnokw(crs) ? nothing : crs
    mappedcrs = isnokw(mappedcrs) ? nothing : mappedcrs
    Projected(parent(l), order, span, sampling, metadata, crs, mappedcrs, dim)
end

GeoInterface.crs(lookup::Projected) = lookup.crs
mappedcrs(lookup::Projected) = lookup.mappedcrs
dim(lookup::Projected) = lookup.dim

@inline function LA.selectindices(l::Projected, sel::LA.Selector; kw...)
    selval = reproject(mappedcrs(l), crs(l), dim(l), val(sel))
    LA._selectindices(l, rebuild(sel; val=selval); kw...)
end
@inline function LA.selectindices(l::Projected, sel::LA.Selector{<:AbstractVector}; kw...)
    selval = reproject(mappedcrs(l), crs(l), dim(l), val(sel))
    LA._selectvec(l, rebuild(sel; val=selval); kw...)
end
@inline function LA.selectindices(l::Projected, sel::LA.IntSelector{<:Tuple}; kw...)
    selval = reproject(mappedcrs(l), crs(l), dim(l), val(sel))
    LA._selecttuple(l, rebuild(sel; val=selval); kw...)
end
@inline LA.selectindices(l::Projected{<:Tuple}, sel::LA.IntSelector{<:Tuple}; kw...) = LA._selectindices(l, sel; kw...)
@inline LA.selectindices(l::Projected{<:Tuple}, sel::LA.IntSelector{<:Tuple{<:Tuple,<:Tuple}}; kw...) = 
    LA._selecttuple(l, sel; kw...)

function LA.selectindices(l::Projected, sel::Between{<:Tuple})
    selval = map(v -> reproject(mappedcrs(l), crs(l), dim(l), v), val(sel))
    LA.between(l, rebuild(sel; val=selval))
end
function LA.selectindices(l::Projected, sel::T) where T<:DD.IntervalSets.Interval
    left, right = map(v -> reproject(mappedcrs(l), crs(l), dim(l), v), (sel.left, sel.right))
    LA.between(l, T(left, right))
end

"""
    Mapped <: AbstractProjected

    Mapped(order, span, sampling, crs, mappedcrs)
    Mapped(; order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(), crs=nothing, mappedcrs)

An [`AbstractSampled`]($DDabssampleddocs) `Lookup`, where the dimension index has
been mapped to another projection, usually lat/lon or `EPSG(4326)`.
`Mapped` matches the dimension format commonly used in netcdf files.

Fields and behaviours are identical to [`Sampled`]($DDsampleddocs) with the addition of
`crs` and `mappedcrs` fields.

The mapped dimension index will be used as for [`Sampled`]($DDsampleddocs),
but to save in another format the underlying `crs` may be used to convert it.
"""
struct Mapped{T,A<:AbstractVector{T},O<:Order,Sp<:Span,Sa<:Sampling,MD,PC<:Union{GeoFormat,Nothing},MC<:Union{GeoFormat,Nothing},D<:Union{AutoDim,Dimension}} <: AbstractProjected{T,O,Sp,Sa}
    data::A
    order::O
    span::Sp
    sampling::Sa
    metadata::MD
    crs::PC
    mappedcrs::MC
    dim::D
end
function Mapped(data=AutoValues();
    order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(),
    metadata=NoMetadata(),
    crs::Union{GeoFormat,Nothing,NoKW}=nokw,
    mappedcrs::Union{GeoFormat,Nothing,NoKW},
    dim=AutoDim()
)
    crs = crs isa NoKW ? nothing : crs
    mappedcrs = mappedcrs isa NoKW ? nothing : mappedcrs
    Mapped(data, order, span, sampling, metadata, crs, mappedcrs, dim)
end
function Mapped(l::Sampled;
    order=order(l), span=span(l), sampling=sampling(l),
    metadata=metadata(l),
    crs::Union{GeoFormat,Nothing}=nothing,
    mappedcrs::Union{GeoFormat,Nothing},
    dim=AutoDim()
)
    Mapped(parent(l), order, span, sampling, metadata, crs, mappedcrs, dim)
end

GeoInterface.crs(lookup::Mapped) = lookup.crs
mappedcrs(lookup::Mapped) = lookup.mappedcrs
dim(lookup::Mapped) = lookup.dim

"""
    convertlookup(dstlookup::Type{<:Lookup}, x)

Convert the dimension lookup between `Projected` and `Mapped`.
Other dimension lookups pass through unchanged.

This is used to e.g. save a netcdf file to GeoTiff.
"""
convertlookup(T::Type{<:Lookup}, A::AbstractDimArray) =
    rebuild(A; dims=convertlookup(T, dims(A)))
convertlookup(T::Type{<:Lookup}, dims::Tuple) = map(d -> convertlookup(T, d), dims)
convertlookup(T::Type{<:Lookup}, d::Dimension) = rebuild(d, convertlookup(T, lookup(d)))
# Non-projected Lookup lookupss pass through
convertlookup(::Type, lookup::Lookup) = lookup
# AbstractProjected passes through if it's the same as dstlookup
convertlookup(::Type{T1}, lookup::T2) where {T1,T2<:T1} = lookup
# Otherwise AbstractProjected needs ArchGDAL
function convertlookup(::Type{<:Mapped}, l::Projected)
    newindex = reproject(crs(l), mappedcrs(l), dim(l), index(l))
    # We use Explicit mode and make a bounds matrix
    # This way the bounds can be saved correctly to NetCDF
    span = if sampling(l) isa Points
        a, b = newindex[1], newindex[end]
        Irregular(LA.isreverse(l) ? (b, a) : (a, b))
    else
        newbounds = reproject(crs(l), mappedcrs(l), dim(l), Dimensions.dim2boundsmatrix(l))
        Explicit(newbounds)
    end
    return Mapped(newindex;
        order=order(l),
        sampling=sampling(l),
        span,
        metadata=metadata(l),
        crs=crs(l),
        mappedcrs=mappedcrs(l),
        dim=dim(l),
    )
end
function convertlookup(::Type{<:Projected}, l::Mapped)
    newindex = _projectedrange(l)
    return Projected(newindex;
        order=order(l),
        span=Regular(step(newindex)),
        sampling=sampling(l),
        metadata=metadata(l),
        crs=crs(l),
        mappedcrs=mappedcrs(l),
        dim=dim(l),
    )
end

_projectedrange(l::Projected) = LinRange(first(l), last(l), length(l))
_projectedrange(l::Mapped) = _projectedrange(span(l), crs(l), l)
function _projectedrange(span, crs, l::Mapped)
    start, stop = reproject(mappedcrs(l), crs, dim(l), [first(l), last(l)])
    LinRange(start, stop, length(l))
end
_projectedrange(::Regular, crs::Nothing, l::Mapped) = LinRange(first(l), last(l), length(l))
function _projectedrange(::T, crs::Nothing, l::Mapped) where T<:Union{Irregular,Explicit}
    error("Cannot convert a Mapped $T index to Projected when crs is nothing")
end

"""
    mappedbounds(x)

Get the bounds converted to the [`mappedcrs`](@ref) value.

Without ArchGDAL loaded, this is just the regular bounds.
"""
function mappedbounds end

mappedbounds(dims::Tuple) = map(mappedbounds, dims)
mappedbounds(dim::Dimension) = mappedbounds(parent(dim), dim)
mappedbounds(::Lookup, dim) = bounds(dim)
mappedbounds(lookup::Projected, dim) = mappedbounds(mappedcrs(lookup), lookup, dim)
mappedbounds(mappedcrs::Nothing, lookup::Projected, dim) =
    error("No mappedcrs attached to $(name(dim)) dimension")
mappedbounds(mappedcrs::GeoFormat, lookup::Projected, dim) =
    _sort(reproject(crs(lookup), mappedcrs, dim, bounds(dim)))

projectedbounds(dims::Tuple) = map(projectedbounds, dims)
projectedbounds(dim::Dimension) = projectedbounds(parent(dim), dim)
projectedbounds(::Lookup, dim) = bounds(dim)
projectedbounds(lookup::Mapped, dim) = projectedbounds(crs(lookup), lookup, dim)
projectedbounds(crs::Nothing, lookup::Mapped, dim) =
    error("No projection crs attached to $(name(dim)) dimension")
projectedbounds(crs::GeoFormat, lookup::Mapped, dim) =
    _sort(reproject(mappedcrs(lookup), crs, dim, bounds(dim)))

_sort((a, b)) = a <= b ? (a, b) : (b, a)

"""
    mappedindex(x)

Get the index value of a dimension converted to the `mappedcrs` value.

Without ArchGDAL loaded, this is just the regular dim value.
"""
function mappedindex end

mappedindex(dims::Tuple) = map(mappedindex, dims)
mappedindex(dim::Dimension) = _mappedindex(parent(dim), dim)

_mappedindex(::Lookup, dim::Dimension) = index(dim)
_mappedindex(lookup::Projected, dim::Dimension) = _mappedindex(mappedcrs(lookup), lookup, dim)
_mappedindex(mappedcrs::Nothing, lookup::Projected, dim) =
    error("No mappedcrs attached to $(name(dim)) dimension")
_mappedindex(mappedcrs::GeoFormat, lookup::Projected, dim) =
    reproject(crs(dim), mappedcrs, dim, index(dim))

projectedindex(dims::Tuple) = map(projectedindex, dims)
projectedindex(dim::Dimension) = _projectedindex(parent(dim), dim)

_projectedindex(::Lookup, dim::Dimension) = index(dim)
_projectedindex(lookup::Mapped, dim::Dimension) = _projectedindex(crs(lookup), lookup, dim)
_projectedindex(crs::Nothing, lookup::Mapped, dim::Dimension) =
    error("No projection crs attached to $(name(dim)) dimension")
_projectedindex(crs::GeoFormat, lookup::Mapped, dim::Dimension) =
    reproject(mappedcrs(dim), crs, dim, index(dim))
