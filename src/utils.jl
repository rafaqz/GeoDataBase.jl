
# File paths, urls and strings

filter_ext(path, ext::AbstractString) =
    filter(fn -> splitext(fn)[2] == ext, readdir(path; join=true))
filter_ext(path, exts::Union{Tuple,AbstractArray}) =
    filter(fn -> splitext(fn)[2] in exts, readdir(path; join=true))
filter_ext(path, ext::Nothing) = readdir(path; join=true)

_maybe_add_suffix(filename::Nothing, suffix) = nothing
_maybe_add_suffix(filename::Nothing, suffix::Union{Nothing,NoKW}) = nothing
_maybe_add_suffix(filename, suffix::Union{Nothing,NoKW}) = filename
function _maybe_add_suffix(filename, suffix)
    base, ext = splitext(filename)
    if string(suffix) == ""
        filename
    else
        return string(base, "_", suffix, ext)
    end
end

# Modified from IsURL.jl, many thanks to @zlatanvasovic
const WINDOWSREGEX = r"^[a-zA-Z]:[\\]"
const URLREGEX = r"^[a-zA-Z][a-zA-Z\d+\-.]*:"

_isurl(str::AbstractString) = !occursin(WINDOWSREGEX, str) && occursin(URLREGEX, str)

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


# Missing values

_missingval_or_missing(x) = _maybe_to_missing(missingval(x))

_maybe_to_missing(::Union{Nothing,NoKW}) = missing
_maybe_to_missing(missingval) = missingval

_writeable_missing(filename::Nothing, T; kw...) = missing
_writeable_missing(filename::AbstractString, T; kw...) = _writeable_missing(T; kw...)
_writeable_missing(::Type{Missing}; verbose=true) = _writeable_missing(UInt8; verbose=true)
function _writeable_missing(T; verbose=true)
    missingval = _type_missingval(T)
    verbose && @info "`missingval` set to $missingval on disk"
    return missingval
end

_type_missingval(::Type{T}) where T = _type_missingval1(Missings.nonmissingtype(T))

_type_missingval1(::Type{T}) where T<:Number = typemin(T)
_type_missingval1(::Type{T}) where T<:Unsigned = typemax(T)
_type_missingval1(::Type{T}) where T<:AbstractString = T("")

_fix_missingval(::Type, ::Union{NoKW,Nothing}) = nothing
_fix_missingval(::AbstractArray, ::Nothing) = nothing
_fix_missingval(A::AbstractArray, ::NoKW) = _fix_missingval(A, Rasters.missingval(A))
_fix_missingval(::AbstractArray{T}, missingval) where T = _fix_missingval(T, missingval)
function _fix_missingval(::Type{T}, missingval::M) where {T,M}
    T1 = nonmissingtype(T)
    if missingval isa T
        missingval
    elseif hasmethod(convert, Tuple{Type{T1},M}) && isreal(missingval) && 
            missingval <= typemax(T1) && missingval >= typemin(T1)
        if T1 <: Integer && !isinteger(missingval) 
            nothing
        else
            convert(T, missingval)
        end
    else
        nothing
    end
end


# Extents

function _extent2dims(to::Extents.Extent; 
    size=nokw, res=nokw, crs=nokw, mappedcrs=nokw, sampling=nokw, kw...
)
    sampling = _match_to_extent(to, isnokw(sampling) ? Intervals(Start()) : sampling)
    _extent2dims(to, size, res; crs, mappedcrs, sampling, kw...)
end
_extent2dims(to::Extents.Extent, size::Union{Nothing,NoKW}, res::Union{Nothing,NoKW}; kw...) =
    throw(ArgumentError("Pass either `size` or `res` keywords or a `Tuple` of `Dimension`s for `to`."))
_extent2dims(to::Extents.Extent, size, res; kw...) = _size_and_res_error()
function _extent2dims(to::Extents.Extent, size::Union{Nothing,NoKW}, res; 
    sampling::Tuple, closed=true, kw...
)
    res = _match_to_extent(to, res)
    ranges = map(values(to), res, sampling) do (start, stop), step, samp
        @assert step >= zero(step) "only positive `res` are supported, got $step"
        if samp isa Intervals
            if locus(samp) isa End
                reverse(range(; start=stop+step, step=-step, stop=start+step))
            else
                r = range(; start, step, stop)
                if locus(samp) isa Start
                    r
                else # Center
                    r .+ step / 2
                end
            end
        else
            range(; start, step, stop)
        end
    end
    return _extent2dims(to, ranges; sampling, kw...)
end
function _extent2dims(to::Extents.Extent, size, res::Union{Nothing,NoKW};
    sampling::Tuple, crs, mappedcrs, closed=true, kw...
)
    size = _match_to_extent(to, size)
    ranges = map(values(to), size, sampling) do (start, stop), length, samp
        r1 = range(; start, stop, length)
        if samp isa Points
            r1
        else
            # We need to buffer extent for a closed interval input so that e.g. 
            # The raster will actually contain points of the extent, as raster
            # pixels are closed/open intervals not closed/closed.
            # Its hard to add a very small amount to any fp number and have 
            # it propagate through to the final extent. But this setup seems to work. 
            # We use the step or r1 to offset the end point of r, the buffer with 10 float
            # steps, which seems to be enough. But it's arbitrary and count be revised.
            nfloatsteps = 10
            step = (stop - start) / length
            if locus(samp) isa End
                start_open = if (closed && start isa AbstractFloat) 
                    prevfloat(start + step, nfloatsteps)
                else
                    start + step
                end
                reverse(range(; start=stop, stop=start_open, length))
            else
                stop_open = if closed && stop isa AbstractFloat 
                    nextfloat(stop - step, nfloatsteps)
                else
                    stop - step
                end
                r = range(; start, stop=stop_open, length)
                if locus(samp) isa Start
                    r
                else # Center
                    r .+ (step / 2)
                end
            end
        end
    end
    return _extent2dims(to, ranges; sampling, crs, mappedcrs)
end
function _extent2dims(::Extents.Extent{K}, ranges; 
    crs, mappedcrs, sampling::Tuple, kw...
) where K
    crs = isnokw(crs) ? nothing : crs 
    mappedcrs = isnokw(mappedcrs) ? nothing : mappedcrs 
    emptydims = map(name2dim, K)
    lookups = map(emptydims, ranges, sampling) do d, range, s
        order = Lookups.orderof(range)
        span = Regular(step(range))
        if d isa SpatialDim && !isnothing(crs)
            Projected(range; sampling=s, order, span, crs, mappedcrs)
        else
            Sampled(range; sampling=s, order, span)
        end
    end
    d = map(rebuild, emptydims, lookups)
    return d
end

function _match_to_extent(::Extents.Extent{K}, x) where K
    if x isa DimTuple 
        map(val, dims(x, map(name2dim, K)))
    elseif x isa NamedTuple 
        values(x[K])
    elseif x isa Tuple 
        x
    else
        map(_ -> x, K)
    end
end

# Geometries

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
        if trait isa GI.FeatureCollectionTrait
            [GI.geometry(f) for f in GI.getfeature(data)]
        elseif isnothing(trait)
            collect(data)
        else
            data
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

# Chunking

# NoKW means true
@inline function _chunks_to_tuple(template, dimorder, chunks::Bool)
    if chunks == true
        if template isa AbstractArray && DA.haschunks(template) == DA.Chunked()
            # Get chunks from the template
            DA.max_chunksize(DA.eachchunk(template))
        else
            # Use defaults
            _chunks_to_tuple(template, dimorder, (X(512), Y(512)))
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

# Constructor helpers

function _raw_check(raw, scaled, missingval, verbose)
    if raw
        # Scaled is false if raw is true
        scaled isa Bool && scaled && verbose && @warn "`scaled=true` set to `false` because of `raw=true`"
        # Only missingval of `nothing` has a meaning with `raw=true`,
        # it turns off missingval completely. Other msissingval values are
        # ignored and a warning is thrown unless verbose=false
        if !isnokwornothing(missingval)
            verbose && @warn "`missingval=$missingval` target value is not used because of `raw=true`"
        end
        return false, isnothing(missingval) ? nothing : Rasters.missingval
    else
        # Otherwise scaled is true and missingval is unchanged
        scaled = isnokw(scaled) ? true : scaled
        return scaled, missingval 
    end
end

function _maybe_drop_single_band(x, dropband::Bool, lazy::Bool)
    dropband || return x
    if hasdim(x, Band()) && size(x, Band()) < 2
         if lazy
             return view(x, Band(1)) # TODO fix dropdims in DiskArrays
         else
             return dropdims(x; dims=Band())
         end
    else
         return x
    end
end


# Memory

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


# Lookups

function _as_intervals(ds::Tuple)
    # Rasterization only makes sense on Sampled Intervals
    interval_dims = map(dims(ds, DEFAULT_POINT_ORDER)) do d
        l = parent(d)
        rebuild(d, rebuild(l; sampling=Intervals(locus(l))))
    end
    return setdims(ds, interval_dims)
end

nolookup_to_sampled(A) = rebuild(A; dims=nolookup_to_sampled(dims(A)))
nolookup_to_sampled(dims::DimTuple) = map(nolookup_to_sampled, dims)
nolookup_to_sampled(d::Dimension) =
    lookup(d) isa NoLookup ? set(d, Sampled(; sampling=Points())) : d


# Metadata

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


# Other

_progress(args...; kw...) = ProgressMeter.Progress(args...; dt=0.1, color=:blue, barlen=50, kw...)

# Function barrier for splatted vector broadcast
@noinline _do_broadcast!(f, x, args...) = broadcast!(f, x, args...)

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

_unwrap(::Val{X}) where X = X
_unwrap(x) = x

# Map filename suffix over a stack
function mapargs(f, st::AbstractRasterStack, args...; kw...)
    layers = map(values(st), args...) do A, mappedargs...
        f(A, mappedargs...; kw...)
    end
    return DD.rebuild_from_arrays(st, Tuple(layers))
end

_without_mapped_crs(f, x) = _without_mapped_crs(f, x, mappedcrs(x))
_without_mapped_crs(f, x, ::Nothing) = f(x)
function _without_mapped_crs(f, dims::DimTuple, mappedcrs::GeoFormat)
    dims1 = setmappedcrs(dims, nothing)
    x = f(dims1)
    if x isa DimTuple
        x = setmappedcrs(x, mappedcrs)
    end
    return x
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


# Warnings and erros

_maybe_warn_replace_missing(replace_missing::NoKW) = nothing
function _maybe_warn_replace_missing(replace_missing)
    @warn """
    `replace_missing` keyword no longer used. Rasters now automatically replaces `missingval` with `missing`. 
    Set `missingval=Rasters.missingval`, to keep the internal missing value, or to replace with some value besides 
    `missing` use e.g. `missingval=Rasters.missingval => NaN` for `NaN`. 
    """
end

@noinline _warn_disk() = @warn "Disk-based objects may be very slow here. User `read` first."

_filenotfound_error(filename) = throw(ArgumentError("file \"$filename\" not found"))

_size_and_res_error() = throw(ArgumentError("Both `size` and `res` keywords are passed, but only one can be used"))

_no_crs_error() = throw(ArgumentError("The provided object does not have a CRS. Use `setcrs` to set one."))


# Rows for extraction etc

# _rowtype returns the complete NamedTuple type for a point row
# This code is entirely for types stability and performance.
# It is used in extract and Rasters.sample
_names(A::AbstractRaster) = (Symbol(name(A)),)
_names(A::AbstractRasterStack) = keys(A)

using DimensionalData.Lookups: _True, _False

_booltype(::_True) = _True()
_booltype(::_False) = _False()
_booltype(x) = x ? _True() : _False()

istrue(::_True) = true
istrue(::_False) = false
istrue(::Type{_True}) = true
istrue(::Type{_False}) = false

# skipinvalid: can G and I be missing. skipmissing: can nametypes be missing
_rowtype(x, g, args...; kw...) = _rowtype(x, typeof(g), args...; kw...)
function _rowtype(
    x, ::Type{G}, i::Type{I} = typeof(size(x)); 
    id=_False(), geometry, index, skipmissing, skipinvalid=skipmissing, names, kw...
) where {G, I}
    _G = istrue(skipinvalid) ? nonmissingtype(G) : G
    _I = istrue(skipinvalid) ? I : Union{Missing, I}
    keys = _rowkeys(id, geometry, index, names)
    types = _rowtypes(x, _G, _I, id, geometry, index, skipmissing, names)
    NamedTuple{keys,types}
end

@generated function _rowtypes(
    x, ::Type{G}, ::Type{I}, id, geometry, index, skipmissing, names::NamedTuple{Names}
) where {G,I,Names}
    ts = Expr(:tuple)
    istrue(id) && push!(ts.args, Int)
    if istrue(skipmissing)
        istrue(geometry) && push!(ts.args, G)
        istrue(index) && push!(ts.args, I)
    else
        istrue(geometry) && push!(ts.args, Union{Missing,G})
        istrue(index) && push!(ts.args, Union{Missing,I})
    end
    :(Tuple{$ts...,_nametypes(x, names, skipmissing)...})
end

@inline _nametypes(::Raster{T}, ::NamedTuple{Names}, sm::_True) where {T,Names} = 
    (nonmissingtype(T),)
@inline _nametypes(::Raster{T}, ::NamedTuple{Names}, sm::_False) where {T,Names} = 
    (Union{Missing,T},)
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

_rowkeys(id::_True, geometry::_False, index::_False, names::NamedTuple{Names}) where Names = (:id, Names...)
_rowkeys(id::_True, geometry::_True, index::_False, names::NamedTuple{Names}) where Names = (:id, :geometry, Names...)
_rowkeys(id::_True, geometry::_True, index::_True, names::NamedTuple{Names}) where Names = (:id, :geometry, :index, Names...)
_rowkeys(id::_True, geometry::_False, index::_True, names::NamedTuple{Names}) where Names = (:id, :index, Names...)
_rowkeys(id::_False, geometry::_False, index::_False, names::NamedTuple{Names}) where Names = Names
_rowkeys(id::_False, geometry::_True, index::_False, names::NamedTuple{Names}) where Names = (:geometry, Names...)
_rowkeys(id::_False, geometry::_True, index::_True, names::NamedTuple{Names}) where Names = (:geometry, :index, Names...)
_rowkeys(id::_False, geometry::_False, index::_True, names::NamedTuple{Names}) where Names = (:index, Names...)