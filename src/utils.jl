
filter_ext(path, ext::AbstractString) =
    filter(fn -> splitext(fn)[2] == ext, readdir(path; join=true))
filter_ext(path, exts::Union{Tuple,AbstractArray}) =
    filter(fn -> splitext(fn)[2] in exts, readdir(path; join=true))
filter_ext(path, ext::Nothing) = readdir(path; join=true)

cleankeys(name) = (_cleankey(name),)
function cleankeys(keys::Union{NamedTuple,Tuple,AbstractArray})
    Tuple(map(_cleankey, keys, ntuple(i -> i, length(keys))))
end

function _cleankey(name::Union{Symbol,AbstractString,Name,NoName}, i=1)
    if name in (NoName(), Symbol(""), Name(Symbol("")))
        Symbol("layer$i")
    else
        Symbol(name)
    end
end

nolookup_to_sampled(A) = rebuild(A; dims=nolookup_to_sampled(dims(A)))
nolookup_to_sampled(dims::DimTuple) = map(nolookup_to_sampled, dims)
nolookup_to_sampled(d::Dimension) =
    lookup(d) isa NoLookup ? set(d, Sampled(; sampling=Points())) : d

function _maybe_use_type_missingval(A::AbstractRaster{T}, source::Source, missingval=nokw) where T
    if ismissing(Rasters.missingval(A))
        newmissingval = missingval isa NoKW ? _type_missingval(Missings.nonmissingtype(T)) : missingval
        A1 = replace_missing(A, newmissingval)
        @warn "`missing` cant be written with $(SOURCE2SYMBOL[source]), missingval for `$(eltype(A1))` of `$newmissingval` used instead"
        return A1
    else
        return A
    end
end

# Create a standardised Metadata object of source T, containing a `Dict{String,Any}`
_metadatadict(s::Source, p1::Pair, pairs::Pair...) = 
    _metadatadict(s, (p1, pairs...))
_metadatadict(::S) where S<:Source = Metadata{S}(Dict{String,Any}())
function _metadatadict(::S, pairs) where S<:Source
    dict = Dict{String,Any}()
    for (k, v) in pairs
        dict[String(k)] = v
    end
    return Metadata{S}(dict)
end

# We often need to convert the locus and the lookup in the same step,
# as doing it in the wrong order can give errors.
# function convert_locus_lookup(M1::Type{<:Lookup}, L1::Type{<:Locus}, dim::Dimension)
#     _convert(S1, L1, sampling(dim), locus(dim), span(dim), dim)
# end

# _convert(::Type{M1}, ::Type{L1}, lookup::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2<:L1} = dim
# _convert(::Type{M1}, ::Type{L1}, lookup::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2} = shiftlocus(L1(), dim)
# _convert(::Type{M1}, ::Type{L1}, lookup::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2} =
#     _convert_by_locus(M1, L1, lookup, l2, span, dim)

# _convert_by_locus(M1, ::Type{Center}, lookup, l2::Union{Start,End}, span, dim) =
#     _convert_by_lookup(M1, dim)
# _convert_by_locus(M1, L1::Type{Union{Start,End}}, lookup, l2::Center, span, dim) =
#     _convert_by_lookup(M1, L1, dim)
# _convert_by_locus(M1, L1::Type{Start}, lookup, l2::End, span, dim) =
#     convertlookup(M1, shiftlocus(L1, dim))
# _convert_by_locus(M1, L1::Type{End}, lookup, l2::Start, span, dim) =
#     convertlookup(M1, shiftlocus(L1, dim))
# _convert_by_locus(M1, ::Type{L1}, lookup, l2::L2, span, dim) where {L1,L2<:L1} =
#     convertlookup(M1, dim)

# # Projected will have an accurate center point on an equal scale, but Mapped may not.
# # So we always shift the locus while in Projected lookup to avoid errors.
# _convert_by_lookup(::Type{Mapped}, dim) = convertlookup(Mapped, shiftlocus(Center(), dim))
# _convert_by_lookup(::Type{Projected}, dim) = shiftlocus(Center(), convertlookup(Projected, dim))


_unwrap(::Val{X}) where X = X
_unwrap(x) = x

_missingval_or_missing(x) = _maybe_nothing_to_missing(missingval(x))

_maybe_nothing_to_missing(::Nothing) = missing
_maybe_nothing_to_missing(missingval) = missingval

maybe_eps(dims::DimTuple; kw...) = map(maybe_eps, dims; kw...)
maybe_eps(dim::Dimension; kw...) = maybe_eps(eltype(dim); kw...)
maybe_eps(x; kw...) = maybe_eps(typeof(x); kw...)
maybe_eps(::Type; kw...) = nothing
maybe_eps(T::Type{<:AbstractFloat}; kw...) = _default_eps(T; kw...)

# These are pretty random defaults, but seem to work
_default_eps(T::Type{<:Float32}; grow=true) = grow ? eps(T) : 100eps(T)
_default_eps(T::Type{<:Float64}; grow=true) = grow ? eps(T) : 1000eps(T)
_default_eps(T::Type{<:Integer}) = T(1)
_default_eps(::Type) = nothing

_writeable_missing(filename::Nothing, T) = missing
_writeable_missing(filename::AbstractString, T) = _writeable_missing(T)
function _writeable_missing(T)
    missingval = _type_missingval(Missings.nonmissingtype(T))
    @info "`missingval` set to $missingval"
    return missingval
end

# Map filename suffix over a stack
function mapargs(f, st::AbstractRasterStack, args...)
    layers = map(values(st), args...) do A, mappedargs...
        f(A, mappedargs...)
    end
    return DD.rebuild_from_arrays(st, Tuple(layers))
end

_without_mapped_crs(f, x) = _without_mapped_crs(f, x, mappedcrs(x))
_without_mapped_crs(f, x, ::Nothing) = f(x)
function _without_mapped_crs(f, dims::DimTuple, mappedcrs::GeoFormat)
    dims1 = setmappedcrs(dims, nothing)
    x = f(dims1)
    return if x isa DimTuple
        setmappedcrs(x, mappedcrs)
    else
        x
    end
end
function _without_mapped_crs(f, A::AbstractRaster, mappedcrs::GeoFormat)
    A = setmappedcrs(A, nothing)
    x = f(A)
    return if x isa AbstractRaster
        setmappedcrs(x, mappedcrs)
    else
        x
    end
end
function _without_mapped_crs(f, st::AbstractRasterStack, mappedcrs::GeoFormat)
    st1 = maplayers(A -> setmappedcrs(A, nothing), st)
    x = f(st1)
    return if x isa AbstractRasterStack
        setmappedcrs(x, mappedcrs)
    else
        x
    end
end

function _extent2dims(to; size=nothing, res=nothing, crs=nothing, kw...) 
    _extent2dims(to, size, res, crs)
end
function _extent2dims(to::Extents.Extent, size::Nothing, res::Nothing, crs)
    isnothing(res) && throw(ArgumentError("Pass either `size` or `res` keywords or a `Tuple` of `Dimension`s for `to`."))
end
function _extent2dims(to::Extents.Extent, size, res, crs)
    isnothing(res) || _size_and_res_error()
end
function _extent2dims(to::Extents.Extent{K}, size::Nothing, res::Real, crs) where K
    tuple_res = ntuple(_ -> res, length(K))
    _extent2dims(to, size, tuple_res, crs)
end
function _extent2dims(to::Extents.Extent{K}, size::Nothing, res, crs) where K
    ranges = map(values(to), res) do bounds, r
        start, stop_closed = bounds
        stop_open = stop_closed + maybe_eps(stop_closed; grow=false)
        length = ceil(Int, (stop_open - start) / r)
        range(; start, step=r, length)
    end
    return _extent2dims(to, ranges, crs)
end
function _extent2dims(to::Extents.Extent{K}, size, res::Nothing, crs) where K
    if size isa Int
        size = ntuple(_ -> size, length(K))
    end
    ranges = map(values(to), size) do bounds, length
        start, stop_closed = bounds
        stop_open = stop_closed + maybe_eps(stop_closed; grow=false)
        step = (stop_open - start) / length
        range(; start, step, length)
    end
    return _extent2dims(to, ranges, crs)
end
function _extent2dims(to::Extents.Extent{K}, ranges, crs) where K
    emptydims = map(name2dim, K)
    lookups = map(ranges) do range
        Projected(range;
            order=ForwardOrdered(),
            sampling=Intervals(Start()),
            span=Regular(step(range)),
            crs,
        )
    end
    d = map(rebuild, emptydims, lookups)
    return d
end

function _as_intervals(ds::Tuple)
    # Rasterization only makes sense on Sampled Intervals
    interval_dims = map(dims(ds, DEFAULT_POINT_ORDER)) do d
        l = parent(d)
        rebuild(d, rebuild(l; sampling=Intervals(locus(l))))
    end
    return setdims(ds, interval_dims)
end

# get geometries from what may be a table with a geometrycolumn or an interable of geometries
# if it has no geometry column and does not iterate valid geometries, error informatively
function _get_geometries(data, ::Nothing)
    # if it's a table, get the geometry column
    geoms = if !(data isa AbstractVector{<:GeoInterface.NamedTuplePoint}) && Tables.istable(data)
        geomcol = first(GI.geometrycolumns(data))
        !in(geomcol, Tables.columnnames(Tables.columns(data))) &&
            throw(ArgumentError("Expected geometries in the column `$geomcol`, but no such column found."))
        isnothing(geomcol) && throw(ArgumentError("No default `geometrycolumn` for this type, please specify it manually."))
        Tables.getcolumn(Tables.columns(data), geomcol)
    elseif data isa AbstractVector
        data
    else
        trait = GI.trait(data)
        if trait isa GI.AbstractFeatureCollectionTrait
            [GI.geometry(f) for f in GI.getfeature(data)]
        elseif trait isa GI.AbstractGeometryCollectionTrait
            GI.getgeom(data)
        elseif trait isa GI.AbstractFeatureTrait
            GI.geometry(data)
        elseif isnothing(trait)
            collect(data)
        elseif trait isa GI.AbstractGeometryTrait
            # data is already a geometry, so return as-is
            data
        else
            ArgumentError("data has $trait, which is not handled")
        end
    end
    # check if data iterates valid geometries before returning
    _check_geometries(geoms)
    return geoms
end
function _get_geometries(data, geometrycolumn::Symbol)
    Tables.istable(data) || throw(ArgumentError("`geometrycolumn` was specified, but `data` is not a table."))
    geoms = Tables.getcolumn(Tables.columns(data), geometrycolumn)
    _check_geometries(geoms)
    return geoms
end
function _get_geometries(data, geometrycolumn::NTuple{<:Any, <:Symbol})
    Tables.istable(data) || throw(ArgumentError("`geometrycolumn` was specified, but `data` is not a table."))
    cols = Tables.columns(data)
    geomcols = (Tables.getcolumn(cols, col) for col in geometrycolumn)
    points = map(geomcols...) do (row...)
        for r in row
            ismissing(r) && return missing
        end
        return row
    end     
    return points
end
function _check_geometries(geoms)
    !isnothing(GI.trait(geoms)) && return
    for g in geoms
        ismissing(g) || !isnothing(GI.geomtrait(g)) || 
            throw(ArgumentError("$g is not a valid GeoInterface.jl geometry"))
    end
    return
end
# to distinguish between objects returned by _get_geometries and other objects
struct IterableOfGeometries end
_warn_disk() = @warn "Disk-based objects may be very slow here. User `read` first."

_filenotfound_error(filename) = throw(ArgumentError("file \"$filename\" not found"))

_progress(args...; kw...) = ProgressMeter.Progress(args...; color=:blue, barlen=50, kw...)

# Function barrier for splatted vector broadcast
@noinline _do_broadcast!(f, x, args...) = broadcast!(f, x, args...)

_size_and_res_error() = throw(ArgumentError("Both `size` and `res` keywords are passed, but only one can be used"))

_no_crs_error() = throw(ArgumentError("The provided object does not have a CRS. Use `setcrs` to set one."))

_type_missingval(::Type{T}) where T = typemin(T)
_type_missingval(::Type{T}) where T<:Unsigned = typemax(T) 

# Modified from IsURL.jl, many thanks to @zlatanvasovic
const WINDOWSREGEX = r"^[a-zA-Z]:[\\]"
const URLREGEX = r"^[a-zA-Z][a-zA-Z\d+\-.]*:"

_isurl(str::AbstractString) = !occursin(WINDOWSREGEX, str) && occursin(URLREGEX, str)

# Run `f` threaded or not, w
function _run(f, range::OrdinalRange, threaded::Bool, progress::Bool, desc::String) 
    p = progress ? _progress(length(range); desc) : nothing
    if threaded
        Threads.@threads :static for i in range
            f(i)
            isnothing(p) || ProgressMeter.next!(p)
        end
    else
        for i in range
            f(i)
            isnothing(p) || ProgressMeter.next!(p)
        end
    end
end

# NoKW means true
@inline function _chunks_to_tuple(template, dims, chunks::Bool)
    if chunks == true
        if template isa AbstractArray && DA.haschunks(template) == DA.Chunked()
            # Get chunks from the template
            DA.max_chunksize(DA.eachchunk(template))
        else
            # Use defaults
            _chunks_to_tuple(template, dims, (X(512), Y(512)))
        end
    else
        nothing
    end
end
@inline function _chunks_to_tuple(template, dimorder, chunks::NTuple{N,Integer}) where N
    n = length(dimorder) 
    if n < N 
        throw(ArgumentError("Length $n tuple needed for `chunks`, got $N"))
    elseif n > N
        (chunks..., ntuple(_ -> 1, Val{n-N}())...)
    else
        chunks
    end
end
@inline function _chunks_to_tuple(template, dimorder, chunks::DimTuple)
    size_one_chunk_axes = map(d -> rebuild(d, 1), otherdims(dimorder, chunks))
    alldims = (chunks..., size_one_chunk_axes...)
    int_chunks = map(val, dims(alldims, dimorder))
    if !isnothing(template)
        if !all(map(>=, size(template), int_chunks))
            @warn "Chunks $int_chunks larger than array size $(size(template)). Using defaults."
            return nothing
        end
    end
    return int_chunks
end
@inline _chunks_to_tuple(template, dimorder, chunks::NamedTuple) =
    _chunks_to_tuple(template, dimorder, DD.kw2dims(chunks))
@inline _chunks_to_tuple(template, dimorder, chunks::Nothing) = nothing
@inline _chunks_to_tuple(template, dims, chunks::NoKW) = nothing


_checkregular(A::AbstractRange) = true
function _checkregular(A::AbstractArray)
    step = stepof(A)
    for i in eachindex(A)[2:end]
        if !(A[i] - A[i-1] ≈ step)
            return false 
        end
    end
    return true
end

function _checkobjmem(obj) 
    f = bytes -> """
        required memory $(bytes) is greater than system memory $(Sys.free_memory()). 
        Use `lazy=true` if you are loading dataset, and only call `read` on a subset after `view`.
        """
    _checkobjmem(f, obj) 
end
_checkobjmem(f, obj) = _checkmem(f, _sizeof(obj))

_checkmem(f, bytes::Int) = Sys.free_memory() > bytes || _no_memory_error(f, bytes)

_sizeof(A::AbstractArray{T}) where T = sizeof(T) * prod(size(A))
_sizeof(st::AbstractRasterStack) = sum(_sizeof, layers(st))
_sizeof(s::AbstractRasterSeries) =
    length(s) == 0 ? 0 : _sizeof(first(s)) * prod(size(s))

function _no_memory_error(f, bytes)
    msg = f(bytes) * """
    If you beleive this is not correct, pass the keyword `checkmem=false` or set `Rasters.checkmem!(false)` 
    and try again. These options may crash your system if the file is actually larger than memory.
    """
    return error(msg)
end


# _rowtype returns the complete NamedTuple type for a point row
# This code is entirely for types stability and performance.
# It is used in extract and Rasters.sample
_names(A::AbstractRaster) = (Symbol(name(A)),)
_names(A::AbstractRasterStack) = keys(A)

using DimensionalData.Lookups: _True, _False
_booltype(x) = x ? _True() : _False()
istrue(::_True) = true
istrue(::_False) = false

# skipinvalid: can G and I be missing. skipmissing: can nametypes be missing
_rowtype(x, g, args...; kw...) = _rowtype(x, typeof(g), args...; kw...)
function _rowtype(
    x, ::Type{G}, i::Type{I} = typeof(size(x)); 
    geometry, index, skipmissing, skipinvalid = skipmissing, names, kw...
) where {G, I}
    _G = istrue(skipinvalid) ? nonmissingtype(G) : G
    _I = istrue(skipinvalid) ? I : Union{Missing, I}
    keys = _rowkeys(geometry, index, names)
    types = _rowtypes(x, _G, _I, geometry, index, skipmissing, names)
    NamedTuple{keys,types}
end


function _rowtypes(
    x, ::Type{G}, ::Type{I}, geometry::_True, index::_True, skipmissing, names::NamedTuple{Names}
) where {G,I,Names}
    Tuple{G,I,_nametypes(x, names, skipmissing)...}
end
function _rowtypes(
    x, ::Type{G}, ::Type{I}, geometry::_True, index::_False, skipmissing, names::NamedTuple{Names}
) where {G,I,Names}
    Tuple{G,_nametypes(x, names, skipmissing)...}
end
function _rowtypes(
    x, ::Type{G}, ::Type{I}, geometry::_False, index::_True, skipmissing, names::NamedTuple{Names}
) where {G,I,Names}
    Tuple{I,_nametypes(x, names, skipmissing)...}
end
function _rowtypes(
    x, ::Type{G}, ::Type{I}, geometry::_False, index::_False, skipmissing, names::NamedTuple{Names}
) where {G,I,Names}
    Tuple{_nametypes(x, names, skipmissing)...}
end

@inline _nametypes(::Raster{T}, ::NamedTuple{Names}, skipmissing::_True) where {T,Names} = (nonmissingtype(T),)
@inline _nametypes(::Raster{T}, ::NamedTuple{Names}, skipmissing::_False) where {T,Names} = (Union{Missing,T},)
# This only compiles away when generated
@generated function _nametypes(
    ::RasterStack{<:Any,T}, ::NamedTuple{PropNames}, skipmissing::_True
) where {T<:NamedTuple{StackNames,Types},PropNames} where {StackNames,Types}
    nt = NamedTuple{StackNames}(map(nonmissingtype, Types.parameters))
    return values(nt[PropNames])
end
@generated function _nametypes(
    ::RasterStack{<:Any,T}, ::NamedTuple{PropNames}, skipmissing::_False
) where {T<:NamedTuple{StackNames,Types},PropNames} where {StackNames,Types}
    nt = NamedTuple{StackNames}(map(T -> Union{Missing,T}, Types.parameters))
    return values(nt[PropNames])
end

_rowkeys(geometry::_False, index::_False, names::NamedTuple{Names}) where Names = Names
_rowkeys(geometry::_True, index::_False, names::NamedTuple{Names}) where Names = (:geometry, Names...)
_rowkeys(geometry::_True, index::_True, names::NamedTuple{Names}) where Names = (:geometry, :index, Names...)
_rowkeys(geometry::_False, index::_True, names::NamedTuple{Names}) where Names = (:index, Names...)
