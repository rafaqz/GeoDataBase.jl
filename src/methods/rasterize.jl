struct _TakeFirst{MV} <: Function
    missingval::MV
end
(tf::_TakeFirst)(a, b) = a === tf.missingval ? b : a
_take_last(a, b) = b

_reduce_op(::typeof(sum)) = Base.add_sum
_reduce_op(::typeof(prod)) = Base.mul_prod
_reduce_op(::typeof(minimum)) = min
_reduce_op(::typeof(maximum)) = max
_reduce_op(::typeof(last)) = _take_last
_reduce_op(f, missingval) = _reduce_op(f)
_reduce_op(::typeof(first), missingval) = _TakeFirst(missingval)
_reduce_op(x) = nothing

_is_op_threadsafe(::typeof(sum)) = true
_is_op_threadsafe(::typeof(prod)) = true
_is_op_threadsafe(::typeof(minimum)) = true
_is_op_threadsafe(::typeof(maximum)) = true
_is_op_threadsafe(f) = false

_reduce_init(reducer, st::AbstractRasterStack, missingval) = map(A -> _reduce_init(reducer, A, missingval), st)
_reduce_init(reducer, ::AbstractRaster{T}, missingval) where T = _reduce_init(reducer, T, missingval)
_reduce_init(reducer, nt::NamedTuple, missingval) = map(x -> _reduce_init(reducer, x, missingval), nt)
_reduce_init(f, x, missingval) = _reduce_init(f, typeof(x), missingval)

_reduce_init(::Nothing, x::Type{T}, missingval) where T = zero(T)
_reduce_init(f::Function, ::Type{T}, missingval) where T = zero(f((zero(nonmissingtype(T)), zero(nonmissingtype(T)))))
_reduce_init(::typeof(sum), ::Type{T}, missingval) where T = zero(nonmissingtype(T))
_reduce_init(::typeof(prod), ::Type{T}, missingval) where T = oneunit(nonmissingtype(T))
_reduce_init(::typeof(minimum), ::Type{T}, missingval) where T = typemax(nonmissingtype(T))
_reduce_init(::typeof(maximum), ::Type{T}, missingval) where T = typemin(nonmissingtype(T))
_reduce_init(::typeof(last), ::Type{T}, missingval) where T = _maybe_to_missing(missingval)

struct FillChooser{F,I,M}
    fill::F
    init::I
    missingval::M
end

struct RasterCreator{E,D,MD,MV,C,MC}
    eltype::E
    to::D
    filename::Union{String,Nothing}
    suffix::String
    name::Symbol
    metadata::MD
    missingval::MV
    crs::C
    mappedcrs::MC
end
function RasterCreator(to::DimTuple; 
    eltype,
    fill,
    missingval,
    filename=nothing,
    suffix="",
    res=nothing, # We shouldn't need this but coverage does
    crs=nothing,
    mappedcrs=nothing,
    name=nothing,
    metadata=Metadata(Dict()),
    kw...
)
    name = Symbol(_filter_name(name, fill))
    to = _as_intervals(to) # Only makes sense to rasterize to intervals
    RasterCreator(eltype, to, filename, suffix, name, metadata, missingval, crs, mappedcrs)
end
RasterCreator(to::AbstractRaster, data; kw...) = RasterCreator(dims(to); kw...)
RasterCreator(to::AbstractRasterStack, data; kw...) = RasterCreator(dims(to); name, kw...)
RasterCreator(to::Nothing, data; kw...) = RasterCreator(_extent(data; kw...); kw...)
RasterCreator(to, data; kw...) = 
    RasterCreator(_extent(to; kw...); kw...)
function RasterCreator(to::Extents.Extent;
    res::Union{Nothing,Real,NTuple{<:Any,<:Real}}=nothing,
    size::Union{Nothing,Int,NTuple{<:Any,Int}}=nothing, 
    crs=nokw,
    mappedcrs=nokw,
    kw...
)
    to_as_dims = _extent2dims(to; size, res, crs, mappedcrs)
    return RasterCreator(to_as_dims; crs, mappedcrs, kw...)
end


# Handles all data validations needed to
# run before rasterizing
struct Rasterizer{T,G,F,R,O,I,M}
    eltype::T
    geom::G
    fillitr::F
    reducer::R
    op::O
    init::I
    missingval::M
    lock::Union{Threads.SpinLock,Nothing}
    shape::Symbol
    boundary::Symbol
    verbose::Bool
    progress::Bool
    threaded::Bool
    threadsafe_op::Bool
end
function Rasterizer(geom, fill, fillitr;
    reducer=nothing,
    op=nothing,
    missingval=nothing,
    shape=nothing,
    eltype=nothing,
    init=nothing,
    boundary=:center,
    filename=nothing,
    verbose=true,
    progress=true,
    threaded=false,
    threadsafe=nothing,
    kw...
)
    # A single geometry does not need a reducing function 
    if !GI.isgeometry(geom)
        isnothing(reducer) && isnothing(op) && !(fill isa Function) && throw(ArgumentError("either reducer, op or fill must be a function"))
    end
 
    op = _reduce_op(reducer)

    threadsafe_op = isnothing(threadsafe) ? _is_op_threadsafe(op) : threadsafe

    shape = if isnothing(shape)
        if GI.isgeometry(geom)
            _geom_shape(geom)
        else
            _geom_shape(first(geom))
        end
    else
        shape
    end

    stable_reductions = (first, last, sum, prod, maximum, minimum)
    if shape == :points &&
        !GI.isgeometry(geom) &&
        !GI.trait(first(geom)) isa GI.PointTrait &&
        !(reducer in stable_reductions)
        @warn "currently `:points` rasterization of multiple non-`PointTrait` geometries may be inaccurate for `reducer` methods besides $stable_reductions. Make a Rasters.jl github issue if you need this to work"
    end
    eltype, missingval, init = get_eltype_missingval(eltype, missingval, fill, fillitr, init, filename, op, reducer)
    lock = threaded ? Threads.SpinLock() : nothing

    return Rasterizer(eltype, geom, fillitr, reducer, op, init, missingval, lock, shape, boundary, verbose, progress, threaded, threadsafe_op)
end
function Rasterizer(data::T; fill, geometrycolumn=nothing, kw...) where T
    Rasterizer(GI.trait(data), data; fill, geometrycolumn, kw...)
end
function Rasterizer(trait::GI.AbstractFeatureCollectionTrait, fc; fill, kw...)
    fillitr = _iterable_fill(trait, fc, fill)
    geoms = map(f -> GI.geometry(f), GI.getfeature(fc))
    Rasterizer(geoms, fill, fillitr; kw...)
end
function Rasterizer(trait::GI.AbstractFeatureTrait, feature; fill, kw...)
    fillitr = _iterable_fill(trait, feature, fill)
    geom = GI.geometry(feature)
    Rasterizer(geom, fill, fillitr; kw...)
end
function Rasterizer(trait::GI.GeometryCollectionTrait, collection; kw...)
    geoms = collect(GI.getgeom(collection))
    Rasterizer(geoms; kw...)
end
function Rasterizer(trait::Nothing, data; fill, geometrycolumn, kw...)
    geoms = _get_geometries(data, geometrycolumn)
    fillitr = _iterable_fill(trait, data, fill)
    Rasterizer(geoms, fill, fillitr; kw...)
end
function Rasterizer(trait::GI.AbstractGeometryTrait, geom; fill, kw...)
    fillitr = _iterable_fill(trait, geom, fill)
    Rasterizer(geom, fill, fillitr; kw...)
end

function get_eltype_missingval(eltype, missingval, fill, fillitr, init, filename, op, reducer)
    filleltype = if fillitr isa NamedTuple
        if all(map(x -> x isa Number, fillitr))
            map(typeof, fillitr)
        elseif all(map(x -> Base.IteratorEltype(x) isa Base.HasEltype, fillitr))
            map(Base.eltype, fillitr)
        else
            map(typeof ∘ first, fillitr) # This is not really correct.
        end
    elseif fill isa Function
        isnothing(init) ? typeof(missingval) : typeof(fill(init))
    elseif Base.IteratorEltype(fillitr) isa Base.HasEltype
        Base.eltype(fillitr)
    else
        typeof(first(fillitr)) # This is not really correct
    end
    init = isnothing(init) ? _reduce_init(reducer, filleltype, missingval) : init
    _get_eltype_missingval(eltype, missingval, filleltype, fillitr, init, filename, op, reducer)
end
function _get_eltype_missingval(eltype, missingval, filleltype, fillitr, init::NamedTuple, filename, op, reducer)
    eltype = eltype isa NamedTuple ? eltype : map(_ -> eltype, init)
    missingval = missingval isa NamedTuple ? missingval : map(_ -> missingval, init)
    filleltype = filleltype isa NamedTuple ? filleltype : map(_ -> filleltype, init)
    em = map(eltype, missingval, filleltype, fillitr, init) do et, mv, fe, fi, i
        _get_eltype_missingval(et, mv, fe, fi, i, filename, op, reducer)
    end
    eltype = map(first, em)
    missingval = map(x -> x[2], em)
    return eltype, missingval, init
end
function _get_eltype_missingval(known_eltype, missingval, filleltype, fillitr, init, filename, op, reducer)
    fillzero = zero(filleltype)
    eltype = if isnothing(known_eltype)
        if fillitr isa Function
            promote_type(typeof(fillitr(init)), typeof(fillitr(fillzero)))
        elseif op isa Function
            promote_type(typeof(op(init, init)), typeof(op(init, fillzero)), typeof(op(fillzero, fillzero)))
        elseif reducer isa Function
            promote_type(typeof(reducer((init, init))), typeof(reducer((init, fillzero))), typeof(reducer((fillzero, fillzero))))
        else
            promote_type(filleltype, typeof(init))
        end
    else
        known_eltype
    end
    missingval = isnothing(missingval) ? _writeable_missing(filename, eltype) : missingval
    # eltype was not the actually array eltype, so promote it with the missingval
    eltype = isnothing(known_eltype) ? promote_type(typeof(missingval), eltype) : eltype
    return eltype, missingval, init
end

_fill_key_error(names, fill) = throw(ArgumentError("fill key $fill not found in table, use one of: $(Tuple(names))"))

# _featurefillval
# Get fill value from a feature, or use fill itself
_featurefillval(feature, fill::Nothing) = first(GI.properties(feature))
_featurefillval(feature, fill::Symbol) = GI.properties(feature)[fill]
_featurefillval(feature, fill::Val) = _featurefillval(feature, _unwrap(fill))
_featurefillval(feature, fill::NamedTuple) = map(f -> _featurefillval(feature, f), _unwrap(fill))
function _featurefillval(feature, fill::NTuple{<:Any,Symbol})
    map(fill) do key
        getproperty(GI.properties(feature), key)
    end |> NamedTuple{fill}
end
_featurefillval(feature, fill) = fill

_filter_name(name, fill::NamedTuple) = keys(fill)
_filter_name(name::NamedTuple, fill::NamedTuple) = keys(fill)
_filter_name(name::Nothing, fill::Nothing) = nothing
_filter_name(name::DimensionalData.NoName, fill::Union{Symbol,NTuple{<:Any,Symbol}}) = fill
_filter_name(name::Union{NamedTuple,Tuple,Array}, fill::NTuple{<:Any,Symbol}) = fill
function _filter_name(name::Union{NamedTuple,Tuple,Array}, fill::Union{Tuple,Array})
    length(name) == length(fill) || throw(ArgumentError("`name` keyword (possibly from `to` object) does not match length of fill. A fix is to use a `NamedTuple` for `fill`."))
    return name isa NamedTuple ? keys(name) : name
end
function _filter_name(name, fill)
    fill isa Union{Symbol,NTuple{<:Any,Symbol}} ? fill : name
end

# A Tuple of `Symbol` is multiple keys to make a RasterStack
_iterable_fill(trait, data, keys::Tuple{Symbol,Vararg}) =
    NamedTuple{keys}(map(k -> _iterable_fill(trait, data, k), keys))
# A Symbol is a Table or FeatureCollection key, it cant be used as fill itself
function _iterable_fill(trait, data, key::Symbol)
    if GI.isfeature(data)
        return get(() -> throw(ArgumentError("feature has no property `:$key`")), GI.properties(data), key)
    elseif trait isa GI.FeatureCollectionTrait
        return [get(() -> throw(ArgumentError("feature has no property `:$key`")), GI.properties(f), key) for f in GI.getfeature(data)]
    end
    cols = Tables.columns(data)
    # For column tables, get the column now
    names = Tables.columnnames(cols)
    key in names || _fill_key_error(names, key)
    return Tables.getcolumn(cols, key)
end
_iterable_fill(trait, data, fill::Function) = fill
_iterable_fill(trait, data, fill::NamedTuple) = begin
    map(f -> _iterable_fill(trait, data, f), fill)
end
function _iterable_fill(trait, data, fill)
    # trait isa Union{GI.AbstractGeometryTrait,GI.FeatureTrait} && return fill
    if trait isa GI.AbstractGeometryTrait || trait isa GI.FeatureTrait
        return fill
    elseif fill isa Number 
        return Iterators.cycle(fill)
    end

    if trait isa GI.FeatureCollectionTrait
        n = GI.nfeature(data)
    elseif Tables.istable(data)
        n = length(Tables.rows(data))
    elseif Base.IteratorSize(data) isa Union{Base.HasShape,Base.HasLength}
        n = length(data)
    else
        return fill
    end
    fillvec = collect(fill)
    l = length(fillvec)
    if l == 1
        # Cycle all length one iterables to fill every row
        return Iterators.cycle(fillvec[1])
    elseif !(l == n)
        throw(ArgumentError("Length of fill $l does not match length of iterator $n"))
    else
        return fillvec
    end
end

_getfill(itrs::NamedTuple, i::Int) = map(itr -> _getfill(itr, i), itrs)
_getfill(itr::AbstractArray, i::Int) = itr[i]
_getfill(itr::Iterators.Cycle, i::Int) = first(itr)
_getfill(itr, i) = itr


const RASTERIZE_KEYWORDS = """
- `fill`: the value or values to fill a polygon with. A `Symbol` or tuple of `Symbol` will
    be used to retrieve properties from features or column values from table rows. An array
    or other iterable will be used for each geometry, in order. `fill` can also be a function of 
    the current value, e.g. `x -> x + 1`.
- `op`: A reducing function that accepts two values and returns one, like `min` to `minimum`.
    For common methods this will be assigned for you, or is not required. But you can use it
    instead of a `reducer` as it will usually be faster.
$GEOM_KEYWORDS
$GEOMETRYCOLUMN_KEYWORD
$PROGRESS_KEYWORD
$VERBOSE_KEYWORD
$THREADED_KEYWORD
- `threadsafe`: specify that custom `reducer` and/or `op` functions are thread-safe, 
    in that the order of operation or blocking does not matter. For example, 
    `sum` and `maximum` are thread-safe, because the answer is approximately (besides
    floating point error) the same after running on nested blocks, or on all the data.
    In contrast, `median` or `last` are not, because the blocking (`median`) or order (`last`) 
    matters.
"""

const RASTERIZE_ARGUMENTS = """
- `reducer`: a reducing function to reduce the fill value for all geometries that
    cover or touch a pixel down to a single value. The default is `last`.
    Any  that takes an iterable and returns a single value will work, including
    custom functions. However, there are optimisations for built-in methods
    including `sum`, `first`, `last`, `minimum`, `maximum`, `extrema` and `Statistics.mean`.
    These may be an order of magnitude or more faster than
    `count` is a special-cased as it does not need a fill value.
$DATA_ARGUMENT
"""

"""
    rasterize([reducer], data; geometrycolumn, kw...)

Rasterize a GeoInterface.jl compatable geometry or feature,
or a Tables.jl table with a `:geometry` column of GeoInterface.jl objects,
or points columns specified by `geometrycolumn`

# Arguments

$RASTERIZE_ARGUMENTS

# Keywords

These are detected automatically from `data` where possible.

$RASTERIZE_KEYWORDS
$FILENAME_KEYWORD
$SUFFIX_KEYWORD

Note on threading. Performance may be much better with `threaded=false`
if `reducer`/`op` are not `threadsafe`. `sum`, `prod`, `maximum`, `minimum`
`count` and `mean` (by combining `sum` and `count`) are threadsafe. If you know
your algorithm is threadsafe, use `threadsafe=true` to allow all optimisations.
Functions passed to `fill` are always threadsafe, and ignore the `threadsafe` argument.

# Example

Rasterize a shapefile for China and plot, with a border.

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates, Shapefile, Downloads
using Rasters.Lookups

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "country_borders.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Load the shapes for china
china_border = Shapefile.Handle(shapefile_name).shapes[10]

# Rasterize the border polygon
china = rasterize(last, china_border; res=0.1, missingval=0, fill=1, boundary=:touches, progress=false)

# And plot
p = plot(china; color=:spring, legend=false)
plot!(p, china_border; fillalpha=0, linewidth=0.6)

savefig("build/china_rasterized.png"); nothing

# output

```

![rasterize](china_rasterized.png)

$EXPERIMENTAL
"""
function rasterize end
function rasterize(reducer::Function, data; kw...)
    rasterize(data; reducer, name=Symbol(string(reducer)), kw...)
end

_count_init_info(init) = @info "`rasterize` with `count` does not use the `init` keyword, $init ignored"
_count_fill_info(fill) = @info "`rasterize` with `count` does not use the `fill` keyword, $fill ignored"
_count_fill(x) = x + 1

# Catch some functions early

# count is faster with an incrementing function as `fill`
function rasterize(reducer::typeof(count), data; fill=nothing, init=nothing, kw...)
    isnothing(init) || _count_init_info(init)
    isnothing(fill) || _count_fill_info(fill)
    rasterize(data; kw..., name=:count, init=0, reducer=nothing, fill=_count_fill, missingval=0)
end
# `mean` is sum ./ count. This is actually optimal with threading, 
# as its means order is irrelevant so its threadsafe.
function rasterize(reducer::typeof(Statistics.mean), data; fill, kw...)
    sums = rasterize(sum, data; kw..., fill)
    counts = rasterize(count, data; kw..., fill=nothing)
    rebuild(sums ./ counts; name=:mean)
end
function rasterize(data; to=nothing, fill, threaded=false, geometrycolumn=nothing, kw...)
    r = Rasterizer(data; fill, threaded, geometrycolumn, kw...)
    rc = RasterCreator(to, data; geometrycolumn, kw..., eltype=r.eltype, fill, missingval=r.missingval)
    allocs = r.shape == :points ? nothing : _burning_allocs(rc.to; threaded)
    return create_rasterize_dest(rc) do dest
        _rasterize!(dest, r; allocs)
    end
end

######################################
# Create a dest Array to rasterize! into

# create_rasterize_dest
# We create a Raster or RasterStack and apply f to it.
# This may be on disk, which is the reason for applying f rather than just
# returning the initialised object - we may need to open it to be able to write.
create_rasterize_dest(f, r::RasterCreator) = create_rasterize_dest(f, r.eltype, r)
# function _create_rasterize_dest(f, dims; fill, name=nothing, init=nothing, kw...)
    # _create_rasterize_dest(f, fill, init, name, dims; fill, kw...)
# end
function create_rasterize_dest(f::Base.Callable, ::NamedTuple{K}, r::RasterCreator) where K
    layers = map(NamedTuple{K}(K), r.missingval, r.eltype) do name, missingval, eltype
        alloc_rasterize(r; eltype, name, missingval, suffix=name) do a
            # We should run `f` here, but it doesn't work yet.
            a
        end
    end
    # Combine layers into a RasterStack
    st = RasterStack(layers)
    # Apply f to the stack. (remove when we can do this in `_alloc_rasterize` while it is open)
    open(f, st; write=true)
    # Return the updated stack
    return st
end
function create_rasterize_dest(f, _, r::RasterCreator)
    result = alloc_rasterize(r) do a
        f(a)
    end
    return result
end

function alloc_rasterize(f, r::RasterCreator;
    eltype=r.eltype,
    name=r.name,
    missingval=r.missingval,
    metadata=r.metadata,
    suffix=r.suffix,
)
    if prod(size(r.to)) == 0  
        throw(ArgumentError("Destination array is is empty, with size $(size(r.to))). Rasterization is not possible"))
    end
    A = create(r.filename, fill=missingval, eltype, r.to; name, missingval, metadata, suffix) do O
        f(O)
    end
    return A
end

"""
    rasterize!([reducer], dest, data; kw...)

Rasterize the geometries in `data` into the [`Raster`](@ref) or [`RasterStack`](@ref) `dest`,
using the values specified by `fill`.

# Arguments

- `dest`: a `Raster` or `RasterStack` to rasterize into.
$RASTERIZE_ARGUMENTS

# Keywords

These are detected automatically from `A` and `data` where possible.

$RASTERIZE_KEYWORDS
$GEOM_KEYWORDS

# Example

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates, Shapefile, GeoInterface, Downloads
using Rasters.Lookups

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "country_borders.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Load the shapes for indonesia
indonesia_border = Shapefile.Handle(shapefile_name).shapes[1]

# Make an empty EPSG 4326 projected Raster of the area of Indonesia
dimz = X(Projected(90.0:0.1:145; sampling=Intervals(), crs=EPSG(4326))),
       Y(Projected(-15.0:0.1:10.9; sampling=Intervals(), crs=EPSG(4326)))

A = zeros(UInt32, dimz; missingval=UInt32(0))

# Rasterize each indonesian island with a different number. The islands are
# rings of a multi-polygon, so we use `GI.getring` to get them all separately.
islands = collect(GeoInterface.getring(indonesia_border))
rasterize!(last, A, islands; fill=1:length(islands), progress=false)

# And plot
p = plot(Rasters.trim(A); color=:spring)
plot!(p, indonesia_border; fillalpha=0, linewidth=0.7)

savefig("build/indonesia_rasterized.png"); nothing

# output

```

![rasterize](indonesia_rasterized.png)

$EXPERIMENTAL
"""
function rasterize! end
rasterize!(reducer::Function, x::RasterStackOrArray, data; kw...) =
    rasterize!(x::RasterStackOrArray, data; reducer, kw...)
function rasterize!(reducer::typeof(count), x::RasterStackOrArray, data; fill=nothing, init=nothing, kw...)
    isnothing(fill) || @info _count_fill_info(fill)
    isnothing(init) || @info _count_init_info(init)
    rasterize!(x::RasterStackOrArray, data; kw..., reducer=nothing, op=nothing, fill=_count_fill, init=0)
end
function rasterize!(x::RasterStackOrArray, data; threaded=false, geometrycolumn=nothing,kw...)
    if prod(size(x)) == 0  
        @warn "Destination is empty, rasterization skipped"
        return x
    end
    r = Rasterizer(data; eltype=eltype(x), threaded, geometrycolumn, kw...)
    allocs = r.shape == :points ? nothing : _burning_allocs(dims(x); threaded)
    return _rasterize!(x, r; allocs)
end
function _rasterize!(A::RasterStackOrArray, r::Rasterizer; allocs=nothing)
    A1 = _prepare_for_burning(A)
    if r.shape == points
        _rasterize_points!(A1, r; allocs)
    else
        _rasterize!(A1, GI.trait(r.geom), r.geom, r.fillitr, r; allocs)
    end
    return A
end
# Single geometry to rasterize
function _rasterize!(A, ::GI.AbstractGeometryTrait, geom, fill, r::Rasterizer; allocs=nothing)
    (; op, init, missingval, lock, shape, boundary, verbose, progress) = r
    if r.shape === :point
        hasburned = _rasterize_points!(A, GI.trait(geom), geom, fill, r)
    else
        ext = _extent(geom)
        V = view(A, Touches(ext))
        # TODO use length when this is fixed in DimensionalData
        prod(size(V)) > 0 || return false

        bools = _init_bools(commondims(V, DEFAULT_POINT_ORDER), BitArray; metadata=metadata(A))
        boolmask!(bools, geom; allocs, lock, shape, boundary, verbose, progress)
        hasburned = any(bools)
        if hasburned
            # Avoid race conditions
            isnothing(lock) || Base.lock(lock)
            _fill!(V, bools, fill, op, init, missingval)
            isnothing(lock) || Base.unlock(lock)
        end
    end
    return hasburned
end
# Fill points
function _rasterize!(A, trait::GI.AbstractPointTrait, point, fill, r::Rasterizer; allocs=nothing)
    return _fill_point!(A, trait, point; fill, lock=r.lock)
end
function _rasterize!(A, trait::Nothing, geoms, fill, r::Rasterizer; allocs=nothing)
    if r.shape === :point
        return _rasterize_points!(A, geoms, fill, r)
    else
        (; reducer, op, fillitr) = r
        # Everything else is rasterized as line or polygon geometries
        return _rasterize_iterable!(A, geoms, reducer, op, fillitr, r, allocs)
    end
end

# We rasterize all iterables from here
function _rasterize_iterable!(A, geoms, reducer, op, fillitr, r::Rasterizer, allocs)
    # Rasterise into a bitarray and reduce it, on one or many threads. 
    # Memory-intensive for large workloads, but thread-safe and safe for e.g. `median`
    if !(fillitr isa Function) && ((r.threaded && !r.threadsafe_op) || isnothing(op))
        (r.threaded && !r.threadsafe_op) && r.verbose && @warn "if `op` is not threadsafe, `threaded=true` may be slower than `threaded=false`"
        return _reduce_bitarray!(reducer, A, geoms, fillitr, r, allocs)
    end
    # Reduce by rasterizing directly on one or many threads,, with a lock
    # TODO: use separate arrays and combine when `r.threadsafe_op == true` ?
    range = _geomindices(geoms)
    burnchecks = _alloc_burnchecks(range)
    _run(range, r.threaded, r.progress, "Rasterizing...") do i
        geom = geoms[i]
        ismissing(geom) && return nothing
        a = _get_alloc(allocs)
        fill = _getfill(fillitr, i)
        burnchecks[i] = _rasterize!(A, GI.trait(geom), geom, fill, r; allocs=a)
        return nothing
    end
    _set_burnchecks(burnchecks, metadata(A), r.verbose)
    return any(burnchecks)
end

################################
# Fast point rasterization
#
# geoms is a iterator of points
_rasterize_points!(A, r::Rasterizer) = _rasterize_points!(A, r.geom, r.fillitr, r)
_rasterize_points!(A, geom, fillitr, r::Rasterizer) =
    _rasterize_points!(A, GI.trait(geom), geom, fillitr, r)

function _rasterize_points!(A, trait::GI.AbstractGeometryTrait, geom, fill, r::Rasterizer)
    points = GI.getpoint(geom)
    fill1 =_iterable_fill(nothing, points, fill)
    _rasterize_points!(A, nothing, points, fill1, r)
end
function _rasterize_points!(A, ::GI.GeometryCollectionTrait, collection, fill, r::Rasterizer)
    # TODO How to handle fill when there is another level of nesting
    hasburned = false
    for geom in _getgeom(collection)
        hasburned |= _rasterize_points!(A, geom, fill, r)
    end
    return hasburned
end
function _rasterize_points!(A, ::Nothing, geoms, fillitr, r::Rasterizer)
    (; reducer, op, missingval, init) = r
    t1 = GI.trait(first(skipmissing(geoms)))
    hasburned = false
    if t1 isa GI.PointTrait
        # Get extent information to properly shift the points
        # to the region of the array during rounding
        ext = Extents.extent(A)
        xrange = ext.X[2] - ext.X[1]
        yrange = ext.Y[2] - ext.Y[1]
        xsize = size(A, X)
        ysize = size(A, Y)
        s = (; ext, xrange, yrange, xsize, ysize)
        return _rasterize_points_inner!(A, geoms, fillitr, s, reducer, op, missingval, init)
    else
        # Recurse down until we hit points
        if r.fillitr isa Function
            for geom in _getgeom(geoms)
                hasburned |= _rasterize_points!(A, GI.trait(geom), geom, fillitr, r)
            end
        else
            if fillitr isa NamedTuple
                ntfill = _maybe_namedtuple_itr(fillitr)
                for (geom, f) in zip(_getgeom(geoms), ntfill)
                    hasburned |= _rasterize_points!(A, GI.trait(geom), geom, f, r)
                end
            else
                for (geom, f) in zip(_getgeom(geoms), fillitr)
                    hasburned |= _rasterize_points!(A, GI.trait(geom), geom, f, r)
                end
            end
        end
        return hasburned
    end
end

@noinline function _rasterize_points_inner!(A, geoms, fillitr::F, s, reducer::R, op::O, missingval, init)::Bool where {F,O,R}
    function xy(p) 
        # TODO handle reversed lookups
        x = round(Int, (GI.x(p) - s.ext.X[1]) / s.xrange * s.xsize) + 1
        y = round(Int, (GI.y(p) - s.ext.Y[1]) / s.yrange * s.ysize) + 1
        (x, y)
    end

    op = reducer == last ? _take_last : op
    if fillitr isa Function
        # We don't need to iterate fill
        points = Iterators.map(xy, geoms)
        return rasterize_points_fillfunc!(fillitr, A, points, missingval, init)
    end

    if op isa Function
        # `Iterators.map` allocates less
        points_fill = Iterators.map(geoms, _maybe_namedtuple_itr(fillitr)) do p, f
            (xy(p), f)
        end
        return rasterize_points_op!(op, A, points_fill, missingval, init)
    else
        # We need to use regular `map` to a vector for sorting later
        points_fill = map(geoms, _maybe_namedtuple_itr(fillitr)) do p, f
            (xy(p), f)
        end
        return rasterize_points_reduce!(reducer, A, points_fill, missingval, init)
    end
end

# Some algorithms don't need sort, like sum
@noinline function rasterize_points_fillfunc!(fillfunc::F, A, points, missingval, init)::Bool where F<:Function
    hasburned = false
    n = 0
    for point in points
        I = dims2indices(A, (X(point[1]), Y(point[2])))
        checkbounds(Bool, A, I...) || continue
         _fill_func!(fillfunc, A, I)
        # Mark that we have written at least one index
        hasburned = true
    end
    return hasburned
end

function _fill_func!(fillfunc, A::Raster, I)
    @inbounds A[I...] = fillfunc(A[I...])
end
function _fill_func!(fillfunc, A::RasterStack, I)
    @inbounds a = A[I...]
    f1 = map(a) do x
        fillfunc(x)
    end
    @inbounds A[I...] = f1
end

# Some reductions don't need sort, like sum
@noinline function rasterize_points_op!(op, A, points_fill, missingval, init)::Bool
    hasburned = false
    n = 0
    
    for (point, fill) in points_fill
        I = dims2indices(A, (X(point[1]), Y(point[2])))
        checkbounds(Bool, A, I...) || continue
         _fill_op!(op, A, fill, init, missingval, I)
        # Mark that we have written at least one index
        hasburned = true
    end
    return hasburned
end

function _fill_op!(op::O, A::Raster, fill, init, missingval, I) where {O<:Function}
    @inbounds a = A[I...]
    f1 = _choose_fill(op, a, FillChooser(fill, init, missingval))
    @inbounds A[I...] = f1
end
function _fill_op!(op::O, A::RasterStack, fill, init, missingval, I) where {O<:Function}
    @inbounds a = A[I...]
    choosers = map(FillChooser, fill, init, missingval)
    f1 = map(a, choosers) do an, fc
        _choose_fill(op, an, fc)
    end
    @inbounds A[I...] = f1
end

type_length(tup::Type{T}) where {T<:Union{Tuple,NamedTuple}} = length(tup.types)

@noinline function rasterize_points_reduce!(reducer, A, points_fill, missingval, init)
    hasburned = false
    # Convert all points to Int
    sort!(points_fill; by=first, alg=Base.Sort.DEFAULT_STABLE)
    prevpoint = first(first(points_fill))
    startind = 1
    for n in eachindex(points_fill)
        point, fill = points_fill[n]
        if prevpoint === point
            continue # We will reduce all of these points together later on
        else
            I = dims2indices(A, (X(prevpoint[1]), Y(prevpoint[2])))
            checkbounds(Bool, A, I...) || continue
            startind = _fill_reduce!(reducer, A, I, points_fill, startind, n, missingval)
        end
        prevpoint = point # Update the previous point to the current
        hasburned = true # Mark that we have written at least one index
    end
    # Fill the last points
    I = dims2indices(A, (X(prevpoint[1]), Y(prevpoint[2])))
    n = lastindex(points_fill) + 1
    checkbounds(Bool, A, I...) && _fill_reduce!(reducer, A, I, points_fill, startind, n, missingval)
    return hasburned
end

# This will not be correct for multiple geometries
function _fill_reduce!(reducer, A, I, pf, startind, i, missingval)
    @inbounds a = A[I...]
    v = _get_fill(reducer, pf[1][2], pf, startind:i - 1)
    x = _reduce_existing(reducer, a, v, missingval)
    @inbounds A[I...] = x
    return i
end

_get_fill(reducer, ::NamedTuple{K}, pf, range) where K = NamedTuple{K}(map(k -> reducer(pf[n][2][k] for n in range), K))
_get_fill(reducer, ::Any, pf, range) = reducer(pf[n][2] for n in range)

_reduce_existing(reducer, as::NamedTuple, vs::NamedTuple, missingvals::NamedTuple) =
    map((a, v, m) -> _reduce_existing(reducer, a, v, m), as, vs, missingvals)
function _reduce_existing(reducer, a, v, missingval)
    if ismissing(missingval) || (!ismissing(a) && a == missingval)
        v
    else
        reducer((a, v)) # This fill fail for mean, median etc
    end
end

# _reduce_bitarray!
#
# Mask `geoms` into each slice of a BitArray with the combined
# dimensions of `x` and `geoms`, then apply the supplied reducing
# function to each cell along the `:geometry` dimension. This uses a
# conditional iterator over the fill values and the corresponding
# rasterized Bool - true included the fill value in the iterator, false excludes it.
#
# We get 64 Bool values to a regular `Int` meaning this doesn't scale too
# badly for large tables of geometries. 64k geometries and a 1000 * 1000
# raster needs 1GB of memory just for the `BitArray`.
# TODO combine these theyre nearly the same
function _reduce_bitarray!(f, st::AbstractRasterStack, geoms, fill::NamedTuple, r::Rasterizer, allocs)
    (; lock, shape, boundary, verbose, progress, threaded) = r
    # Define mask dimensions, the same size as the spatial dims of x
    spatialdims = commondims(st, DEFAULT_POINT_ORDER)
    # Mask geoms as separate bool layers
    masks = boolmask(geoms; to=st, collapse=false, metadata=metadata(st), allocs, lock, shape, boundary, verbose, progress)
    # Use a generator over the array axis in case the iterator has no length
    geom_axis = axes(masks, Dim{:geometry}())
    fill = map(itr -> [v for (_, v) in zip(geom_axis, itr)], fill)
    T = eltype(st)
    range = axes(st, Y())
    _run(range, threaded, progress, "Reducing...") do y
        _reduce_bitarray_loop(f, st, T, fill, masks, y)
    end
end
function _reduce_bitarray!(f, A::AbstractRaster, geoms, fill, r::Rasterizer, allocs)
    (; lock, shape, boundary, verbose, progress, threaded) = r
    # Define mask dimensions, the same size as the spatial dims of x
    spatialdims = commondims(A, DEFAULT_POINT_ORDER)
    # Mask geoms as separate bool layers
    masks = boolmask(geoms; to=A, collapse=false, metadata=metadata(A), lock, shape, boundary, verbose, progress, threaded)
    # Use a generator over the array axis in case the iterator has no length
    geom_axis = parent(axes(masks, Dim{:geometry}()))
    fill = [val for (i, val) in zip(geom_axis, fill)]
    T = eltype(A)
    range = axes(A, Y())
    _run(range, threaded, progress, "Reducing...") do y
        _reduce_bitarray_loop(f, A, T, fill, masks, y)
    end
    return A
end

function _reduce_bitarray_loop(f, A, ::Type{T}, fill, masks, y) where T
    for x in axes(A, X())
        D = (X(x), Y(y))
        # Do DimensionalData.jl indexing manually to avoid taking a view of the index and TwicePrecision problems
        I = dims2indices(masks, D)
        newval = _apply_reduction!(T, f, fill, view(parent(masks), I...))::Union{T,Nothing}
        if !isnothing(newval)
            @inbounds A[D...] = newval::T
        end
    end
end

# _apply_reduction!
#
# Apply a reducing function over an iterable
# This is applied for all reducing methods that don't have a matching `op` method
@inline function _apply_reduction!(::Type{T}, f, fill, pixel_geom_list) where T
    any(pixel_geom_list) || return nothing
    iterator = (fl for (fl, b) in zip(fill, pixel_geom_list) if b && !ismissing(fl))
    return convert(T, f(iterator))
end
@inline function _apply_reduction!(::Type{T}, f, fill::NamedTuple, pixel_geom_list) where T
    any(pixel_geom_list) || return nothing
    vals = map(fill) do fill
        iterator = (fl for (fl, b) in zip(fill, pixel_geom_list) if b && !ismissing(fl))
        f(iterator)
    end
    return convert(T, vals)
end


######################################
# Fill

_fill!(st::AbstractRasterStack, B, fill::NamedTuple, op, init::Union{Nothing,NamedTuple}, missingval) =
    _fill!(st, B, fill, op, init, map(_ -> missingval, fill))
function _fill!(st::AbstractRasterStack, B, fill::NamedTuple, op, init::Union{Nothing,NamedTuple}, missingval::NamedTuple)
    init = isnothing(init) ? map(_ -> nothing, fill) : init
    foreach(DimensionalData.layers(st), fill, init, missingval) do a, f, i, mv
        _fill!(a, B, f, op, i, mv)
    end
    return st
end
# If the array is initialised, we can use the existing values
function _fill!(A::AbstractRaster{T}, B, fill, op, init, missingval) where T
    broadcast_dims!(A, A, B) do a, b
        convert(T, b ? _choose_fill(op, a, FillChooser(fill, init, missingval)) : a)::T
    end
    return A
end

# No op fill is a function, apply it unless missingval
Base.@assume_effects :total _choose_fill(op::Nothing, a, fc::FillChooser{<:Function}) =
    a == missingval ? fc.fill(fc.init) : fc.fill(a)
# No op fill is a function, apply it unless missing
Base.@assume_effects :total _choose_fill(op::Nothing, a, fc::FillChooser{<:Function,<:Any,Missing}) =
    ismissing(a) ? fc.fill(fc.init) : fc.fill(a)
# No op fill is a function, no init
Base.@assume_effects :total _choose_fill(op::Nothing, a, fc::FillChooser{<:Function,Nothing}) = fc.fill(a)
# No op fill is a function, no init fill a (repeated to avoid ambiguity)
Base.@assume_effects :total _choose_fill(op::Nothing, a, fc::FillChooser{<:Function,Nothing,Missing}) = fc.fill(a)
# Op is a function, fill is not, missingval===missing
# apply reducing op to a and fill, or to init and fill if a equals missing and init exists
Base.@assume_effects :total function _choose_fill(op::F, a, fc::FillChooser{<:Any,<:Any,Missing}) where F<:Function
    a1 = ismissing(a) ? fc.init : a
    _apply_op(op, a1, fc.fill)
end
Base.@assume_effects :total function _choose_fill(op::F, a, fc::FillChooser) where F<:Function
    a1 = a === fc.missingval ? fc.init : a
    _apply_op(op, a1, fc.fill)
end
Base.@assume_effects :total function _choose_fill(op::F, a, fc::FillChooser{<:Any,Nothing,Missing}) where F<:Function
    _apply_op(op, a, fc.fill)
end
# No op, fill is a value - this is just one geometry
Base.@assume_effects :total function _choose_fill(op, a, fc::FillChooser)
    fc.fill
end
# Op is a function, fill is not, missingval===missing
# apply reducing op to a and fill, or to init and fill if a equals missingval and init exists
# @inline function _choose_fill(a, fill, op::F, init::Nothing, missingval) where F<:Function
#     _apply_op(op, a, fill)
# end

# apply reducing op to current value and fill value
Base.@assume_effects :total _apply_op(op::Nothing, a1, fill) = fill
Base.@assume_effects :total _apply_op(op::F, a1, fill) where F<:Function = op(a1, fill)


_maybe_namedtuple_itr(nt::NamedTuple{K}) where K = 
    (NamedTuple{K}(xs) for xs in zip(nt...)) 
_maybe_namedtuple_itr(itr) = itr
