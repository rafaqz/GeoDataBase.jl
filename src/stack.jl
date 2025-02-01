# Accept either Symbol or String keys, but always convert to Symbol
const Key = Union{Symbol,AbstractString}

const MAX_STACK_SIZE = 200

"""
    AbstractRasterStack

Abstract supertype for objects that hold multiple [`AbstractRaster`](@ref)s
that share spatial dimensions.

They are `NamedTuple`-like structures that may either contain `NamedTuple`
of [`AbstractRaster`](@ref)s, string paths that will load [`AbstractRaster`](@ref)s,
or a single path that points to a file containing multiple layers, like
NetCDF or HDF5. Use and syntax is similar or identical for all cases.

`AbstractRasterStack` can hold layers that share some or all of their dimensions.
They cannot have the same dimension with different length or spatial extent as
another layer.

`getindex` on an `AbstractRasterStack` generally returns a memory backed standard
[`Raster`](@ref). `raster[:somelayer] |> plot` plots the layers array,
while `raster[:somelayer, X(1:100), Band(2)] |> plot` will plot the
subset without loading the whole array.

`getindex` on an `AbstractRasterStack` with a key returns another stack with
`getindex` applied to all the arrays in the stack.
"""
abstract type AbstractRasterStack{K,T,N,L} <: AbstractDimStack{K,T,N,L} end

missingval(stack::AbstractRasterStack) = getfield(stack, :missingval)
missingval(s::AbstractRasterStack, name::Symbol) = _singlemissingval(missingval(s), name)

filename(stack::AbstractRasterStack{<:Any,<:Any,<:Any,<:NamedTuple}) = 
    map(s -> filename(s), layers(stack))
filename(stack::AbstractRasterStack{<:Any,<:Any,<:Any,<:Union{FileStack,OpenStack}}) = 
    filename(parent(stack))

isdisk(st::AbstractRasterStack) = any(isdisk, layers(st))

setcrs(x::AbstractRasterStack, crs) = set(x, setcrs(dims(x), crs)...)
setmappedcrs(x::AbstractRasterStack, mappedcrs) = set(x, setmappedcrs(dims(x), mappedcrs)...)

_singlemissingval(mvs::NamedTuple, name) = mvs[name]
_singlemissingval(mv, name) = mv

function _maybe_collapse_missingval(mvs::NamedTuple)
    mv1, mvs_rest = Iterators.peel(mvs)
    for mv in mvs_rest
        mv === mv1 || return mvs
    end
    return mv1
end
_maybe_collapse_missingval(::NoKW) = nothing
_maybe_collapse_missingval(mv) = mv

# DimensionalData methods ######################################################

# Always read a stack before loading it as a table.
DD.DimTable(stack::AbstractRasterStack; checkmem=CHECKMEM[]) =
    invoke(DD.DimTable, Tuple{DD.AbstractDimStack}, read(stack; checkmem))

function DD.layers(s::AbstractRasterStack{<:Any,<:Any,<:Any,<:FileStack{<:Any,Keys}}) where Keys
    NamedTuple{Keys}(map(K -> s[K], Keys))
end
function DD.layers(s::AbstractRasterStack{<:Any,<:Any,<:Any,<:OpenStack{<:Any,Keys}}) where Keys
    NamedTuple{Keys}(map(K -> s[K], Keys))
end

function DD.rebuild(
    s::AbstractRasterStack,
    data,
    dims=dims(s),
    refdims=refdims(s),
    layerdims=DD.layerdims(s),
    metadata=metadata(s),
    layermetadata=DD.layermetadata(s),
    missingval=missingval(s),
)
    DD.basetypeof(s)(data, dims, refdims, layerdims, metadata, layermetadata, missingval)
end
function DD.rebuild(s::AbstractRasterStack;
    data=parent(s),
    dims=dims(s),
    refdims=refdims(s),
    layerdims=DD.layerdims(s),
    metadata=metadata(s),
    layermetadata=DD.layermetadata(s),
    missingval=missingval(s),
)
    DD.basetypeof(s)(
        data, dims, refdims, layerdims, metadata, layermetadata, missingval
    )
end

function DD.rebuild_from_arrays(
    s::AbstractRasterStack{<:Union{FileStack{<:Any,Keys},OpenStack{<:Any,Keys}}}, 
    das::Tuple{Vararg{AbstractDimArray}}; 
    kw...
) where Keys
    DD.rebuild_from_arrays(s, NamedTuple{Keys}(das); kw...)
end
function DD.rebuild_from_arrays(
    s::AbstractRasterStack, 
    das::NamedTuple{<:Any,<:Tuple{Vararg{AbstractDimArray}}};
    refdims=refdims(s),
    metadata=DD.metadata(s),
    dims=nothing,
    layerdims=map(DD.basedims, das),
    layermetadata=map(DD.metadata, das),
    missingval=_maybe_collapse_missingval(map(missingval, das)),
)
    data = map(parent, das)
    if isnothing(dims)
        # invokelatest avoids compiling this for other paths
        dims = Base.invokelatest() do
            DD.combinedims(collect(das))
        end
        rebuild(s; data, dims, refdims, layerdims, metadata, layermetadata, missingval)
    else
        rebuild(s; data, dims, refdims, layerdims, metadata, layermetadata, missingval)
    end
end

# Base methods #################################################################

DD.name(s::AbstractRasterStack) = keys(s)
Base.names(s::AbstractRasterStack) = keys(s)
Base.copy(stack::AbstractRasterStack) = maplayers(copy, stack)

#### Stack getindex ####
# Different to DimensionalData as we construct a Raster
Base.@constprop :aggressive @propagate_inbounds Base.getindex(s::AbstractRasterStack, name::AbstractString) = 
    s[Symbol(name)]
Base.@constprop :aggressive @propagate_inbounds function Base.getindex(s::AbstractRasterStack, name::Symbol)
    data_ = parent(s)[name]
    dims_ = dims(s, DD.layerdims(s, name))
    metadata = DD.layermetadata(s, name)
    mv = missingval(s, name)
    mv = mv isa eltype(data_) ? mv : nothing
    Raster(data_, dims_, refdims(s), name, metadata, mv)
end

# Concrete AbstrackRasterStack implementation #################################################

"""
    RasterStack <: AbstrackRasterStack

    RasterStack(data...; name, kw...)
    RasterStack(data::Union{Vector,Tuple}; name, kw...)
    RasterStack(data::NamedTuple; kw...))
    RasterStack(data::RasterStack; kw...)
    RasterStack(data::Raster; layersfrom=Band, kw...)
    RasterStack(filepath::AbstractString; kw...)

Load a file path or a `NamedTuple` of paths as a `RasterStack`, or convert arguments, a
`Vector` or `NamedTuple` of `Raster`s to `RasterStack`.

# Arguments

- `data`: A `NamedTuple` of [`Raster`](@ref)s or `String`, or a `Vector`, `Tuple` or splatted
    arguments of [`Raster`](@ref). The latter options must pass a `name` keyword argument.
- `filepath`: A file (such as netcdf or tif) to be loaded as a stack, or a directory path
    containing multiple files.

# Keywords

- `name`: Used as stack layer names when a `Tuple`, `Vector` or splat of `Raster` is passed in.
    Has no effect when `NameTuple` is used - the `NamedTuple` keys are the layer names.
$GROUP_KEYWORD
- `metadata`: A `Dict` or `DimensionalData.Metadata` object.
$MISSINGVAL_KEYWORD
    For `RasterStack` a `NamedTuple` can also be passed if layers
    should have different `missingval`.
$CONSTRUCTOR_CRS_KEYWORD
$CONSTRUCTOR_MAPPEDCRS_KEYWORD
- `refdims`: `Tuple` of `Dimension` that the stack was sliced from.

For when one or multiple filepaths are used:

$DROPBAND_KEYWORD
$LAZY_KEYWORD
$RAW_KEYWORD
$SCALED_KEYWORD
$SOURCE_KEYWORD

For when a single `Raster` is used:

- `layersfrom`: `Dimension` to source stack layers from if the file is not already multi-layered.
    `nothing` is default, so that a single `RasterStack(raster)` is a single layered stack.
    `RasterStack(raster; layersfrom=Band)` will use the bands as layers.

```julia
files = (temp="temp.tif", pressure="pressure.tif", relhum="relhum.tif")
stack = RasterStack(files; mappedcrs=EPSG(4326))
stack[:relhum][Lat(Contains(-37), Lon(Contains(144))
```
"""
struct RasterStack{K,T,N,L<:Union{FileStack,OpenStack,NamedTuple},D<:Tuple,R<:Tuple,LD<:NamedTuple,M,LM,MV} <: AbstractRasterStack{K,T,N,L}
    data::L
    dims::D
    refdims::R
    layerdims::LD
    metadata::M
    layermetadata::LM
    missingval::MV
end
function RasterStack{K,T,N}(
    data::L, dims::D, refdims::R, layerdims::LD, metadata::Me, layermetadata::LM, missingval::Mi
) where {K,T,N,L,D,R,LD<:NamedTuple{K},Me,LM,Mi}
    RasterStack{K,T,N,L,D,R,LD,Me,LM,Mi}(data, dims, refdims, layerdims, metadata, layermetadata, missingval)
end
function RasterStack(
    data, dims, refdims, layerdims::LD, metadata, layermetadata, missingval
) where LD<:NamedTuple{K} where K
    T = DD.data_eltype(data)
    N = length(dims)
    RasterStack{K,T,N}(data, dims, refdims, layerdims, metadata, layermetadata, missingval)
end
function RasterStack(
    data::Union{FileStack,OpenStack,NamedTuple};
    dims::Tuple,
    refdims::Tuple=(),
    layerdims::NamedTuple,
    metadata=nokw,
    layermetadata=nokw,
    missingval=nokw,
    kw...
)
    K = keys(data)
    # Handle values that musbe be `NamedTuple`
    layermetadata = if layermetadata isa NamedTuple
        layermetadata
    elseif layermetadata isa Union{Nothing,NoKW,NoMetadata}
        NamedTuple{K}(map(_ -> NoMetadata(), K))
    else
        throw(ArgumentError("$layermetadata is not a valid input for `layermetadata`. Try a `NamedTuple` of `Dict`, `MetaData` or `NoMetadata`"))
    end
    metadata = isnokw(metadata) ? NoMetadata() : metadata
    missingval = _maybe_collapse_missingval(missingval)
    layerdims = if layerdims isa NamedTuple 
        layerdims 
    else
        NamedTuple{K}(ntuple(i -> layerdims[i], Val{length(K)}()))
    end
    st = RasterStack(
        data, dims, refdims, layerdims, metadata, layermetadata, missingval
    )
    return _postprocess_stack(st, dims; kw...)
end
# Convert Tuple/Array of array to NamedTuples using name/key
function RasterStack(data::Tuple{Vararg{AbstractArray}}, dims::Tuple;
    name::Union{Tuple,AbstractArray,NamedTuple,Nothing}=nokw,
    kw...
)
    isnokw(name) && throw(ArgumentError("Pass a Tuple, Array or NamedTuple of names to the `name` keyword"))
    return RasterStack(NamedTuple{cleankeys(name)}(data), dims; kw...)
end
# Multi Raster stack from NamedTuple of AbstractArray
function RasterStack(data::NamedTuple{<:Any,<:Tuple{Vararg{AbstractArray}}}, dims::Tuple; 
    layerdims=nokw,
    kw...
)
    if isnokw(layerdims)
        # TODO: make this more sophisticated and match dimension length to axes?
        # We don't worry about Raster keywords because these rasters will be deconstructed
        # again later, and `kw` will define the RasterStack keywords
        layers = map(data) do A
            Raster(A, dims[1:ndims(A)])
        end
        return RasterStack(layers; kw...)
    else
        return RasterStack(data; dims, layerdims, kw...)
    end
end
# Multi Raster stack from AbstractDimArray splat
RasterStack(layers::AbstractDimArray...; kw...) = RasterStack(layers; kw...)
# Multi Raster stack from tuple with `name` keyword
function RasterStack(layers::Tuple{Vararg{AbstractDimArray}};
    name=map(name, layers),
    kw...
)
    RasterStack(NamedTuple{cleankeys(name)}(layers); kw...)
end
# Multi RasterStack from NamedTuple
# This method is called after most other RasterStack methods.
function RasterStack(layers::NamedTuple{K,<:Tuple{Vararg{AbstractDimArray}}};
    resize::Union{Function,NoKW}=nokw,
    _layers=resize isa NoKW ? layers : resize(layers),
    dims::Tuple=DD.combinedims(_layers...),
    refdims::Tuple=(),
    missingval=map(missingval, _layers),
    metadata=NoMetadata(),
    layermetadata::NamedTuple{K}=map(DD.metadata, _layers),
    layerdims::NamedTuple{K}=map(DD.basedims, _layers),
    kw...
) where K
    data = map(parent, _layers)
    st = RasterStack(data;
        dims, refdims, layerdims, metadata, layermetadata, missingval
    )
    return _postprocess_stack(st, dims; kw...)
end
# Stack from table and dims args
function RasterStack(table, dims::Tuple;
    name=_not_a_dimcol(table, dims),
    kw...
)
    Tables.istable(table) || throw(ArgumentError("object $(typeof(table)) is not a valid input to `RasterStack`"))
    if name isa Symbol
        col = Tables.getcolumn(table, name)
        layers = NamedTuple{(name,)}((reshape(col, map(length, dims)),))
    else
        layers = map(name) do k
            col = Tables.getcolumn(table, k)
            reshape(col, map(length, dims))
        end |> NamedTuple{cleankeys(name)}
    end
    return RasterStack(layers, dims; kw...)
end
# Stack from a Raster
function RasterStack(A::AbstractDimArray;
    layersfrom=nokw,
    name=nokw,
    refdims::Tuple=refdims(A),
    metadata=metadata(A),
    missingval=missingval(A),
    kw...
)
    name = name isa Union{AbstractString,Symbol,Name} ? (name,) : name
    layers = if isnokw(layersfrom)
        name = if isnokw(name)
            name = DD.name(A) in (NoName(), Symbol(""), Name(Symbol(""))) ? ("layer1",) : DD.name(A)
        else
            name
        end
        NamedTuple{cleankeys(name)}((A,))
    else
        name = isnokw(name) ? _layerkeysfromdim(A, layersfrom) : name
        slices = slice(A, layersfrom)
        NamedTuple{cleankeys(name)}(slices)
    end
    return RasterStack(layers; refdims, metadata, missingval, kw...)
end
# Stack from stack and dims args
RasterStack(st::DD.AbstractDimStack, dims::Tuple; kw...) = RasterStack(st; dims, kw...)
# RasterStack from another stack
function RasterStack(s::DD.AbstractDimStack;
    data=parent(s),
    dims::Union{Tuple,NoKW}=dims(s),
    refdims::Tuple=refdims(s),
    metadata=metadata(s),
    layerdims=DD.layerdims(s),
    layermetadata=DD.layermetadata(s),
    missingval=missingval(s),
    kw...
)
    return RasterStack(
        data, dims; refdims, layerdims, metadata, layermetadata, missingval, kw...
    )
end
# Multi-file stack from strings
function RasterStack(
    filenames::Union{AbstractArray{<:AbstractString},Tuple{<:AbstractString,Vararg}};
    name=map(filekey, filenames),
    kw...
)
    RasterStack(NamedTuple{cleankeys(name)}(filenames); kw...)
end
function RasterStack(filenames::NamedTuple{K,<:Tuple{<:AbstractString,Vararg}};
    source=nokw,
    metadata=nokw,
    resize=nokw,
    layermetadata::Union{NoKW,NamedTuple{K}}=nokw,
    layerdims::Union{NoKW,NamedTuple{K}}=nokw,
    missingval=nokw,
    replace_missing=nokw,
    scaled=nokw,
    raw=false,
    verbose=true,
    kw...
) where K
    _maybe_warn_replace_missing(replace_missing)
    scaled, missingval = _raw_check(raw, scaled, missingval, verbose)

    layermissingval = collect(_stack_nt(filenames, missingval))
    fn = collect(filenames)
    layermetadata = layermetadata isa NamedTuple ? collect(layermetadata) : map(_ -> NoKW(), fn)
    layerdims = layerdims isa NamedTuple ? collect(layerdims) : map(_ -> NoKW(), fn)
    layers = map(K, fn, layermetadata, layerdims, layermissingval) do name, fn, md, d, mv
        Raster(fn; 
            source=sourcetrait(fn, source), 
            dims=d, name, metadata=md, missingval=mv, scaled, verbose, kw...
       )
    end
    return RasterStack(NamedTuple{K}(layers); resize, metadata)
end
# Stack from a String
function RasterStack(filename::AbstractString;
    lazy::Bool=false,
    dropband::Bool=true,
    raw::Bool=false,
    source::Union{Symbol,Source,NoKW}=nokw,
    missingval=nokw,
    name=nokw,
    group::Union{Symbol,AbstractString,NoKW}=nokw,
    scaled::Union{Bool,NoKW}=nokw,
    coerce=nokw,
    verbose::Bool=true,
    replace_missing=nokw, # deprecated
    kw...
)
    _maybe_warn_replace_missing(replace_missing)
    scaled, missingval = _raw_check(raw, scaled, missingval, verbose)

    source = sourcetrait(filename, source)
    st = if isdir(filename) && !(source isa Zarrsource)
        # Load as a whole directory
        filenames = readdir(filename)
        length(filenames) > 0 || throw(ArgumentError("No files in directory $filename"))
        # Detect keys from names
        name = if isnokw(name)
            stripped = lstrip.(x -> x in (" ", "_"), (x -> x[1:end]).(filenames))
            Symbol.(replace.(first.(splitext.(stripped)), Ref(" " => "_")))
        else
            name
        end
        RasterStack(joinpath.(Ref(filename), filenames);
            missingval, scaled, coerce, lazy, dropband, group, kw...
        )
    else
        # Load as a single file
        if haslayers(source)
            # With multiple named layers
            l_st = _layer_stack(filename; 
                source, name, lazy, group, missingval, scaled, coerce, kw...
            )
            # Maybe split the stack into separate arrays to remove extra dims.
            if isnokw(name)
                l_st
            else
                maplayers(identity, l_st)
            end
        else
            # With bands actings as layers
            raster = Raster(filename; 
                source, lazy, missingval, scaled, coerce, dropband=false,
            )
            RasterStack(raster; kw...)
        end
    end

    # Maybe drop the Band dimension
    return _maybe_drop_single_band(st, dropband, lazy)
end

# TODO test this properly
function DD.modify(f, s::AbstractRasterStack{<:FileStack{<:Any,K}}) where K
    data = open(s) do o
        map(K) do k
            f(parent(ost)[k])
        end
    end
    rebuild(s; data)
end

# Open a single file stack
function Base.open(f::Function, st::AbstractRasterStack{K,T,<:Any,<:FileStack{X}}; kw...) where {X,K,T}
    ost = OpenStack{X,K,T}(parent(st))
    # TODO is this needed?
    layers = map(K) do k
        ost[k]
    end |> NamedTuple{K}
    out = f(rebuild(st; data=layers))
    close(ost)
    return out
end
# Open a multi-file stack or just apply f to a memory backed stack
function Base.open(f::Function, st::AbstractRasterStack{<:Any,<:Any,<:Any,<:NamedTuple}; kw...)
    isdisk(st) ? _open_layers(f, st) : f(st)
end

# Open all layers through nested closures, applying `f` to the rebuilt open stack
_open_layers(f, st) = _open_layers(f, st, DD.layers(st), NamedTuple())
function _open_layers(f, st, unopened::NamedTuple{K}, opened::NamedTuple) where K
    open(first(unopened)) do open_layer
        layer_nt = NamedTuple{(first(K),)}((open_layer,))
        _open_layers(f, st, Base.tail(unopened), merge(opened, layer_nt))
    end
end
function _open_layers(f, st, unopened::NamedTuple{()}, opened::NamedTuple)
    f(rebuild(st; data=opened))
end

function _postprocess_stack(st, dims;
    crs=nokw,
    mappedcrs=nokw,
    name=nokw,
    resize=nokw,
)
    dims = format(dims, CartesianIndices(st))
    dims = isnokw(crs) ? dims : setcrs(dims, crs)
    dims = isnokw(mappedcrs) ? dims : setmappedcrs(dims, mappedcrs)
    st = rebuild(st; dims)
    return isnokw(name) ? st : st[Dimensions._astuple(name)]
end

function _layer_stack(filename;
    source=nokw,
    dims=nokw,
    refdims=(),
    name=nokw,
    group=nokw,
    metadata=nokw,
    layermetadata=nokw,
    layerdims=nokw,
    missingval=nokw,
    crs=nokw,
    mappedcrs=nokw,
    coerce=convert,
    scaled=nokw,
    checkmem=true,
    lazy=false,
    kw...
)
    data, field_kw = _open(filename; source) do ds
        layers = _layers(ds, name, group)
        # Create a Dict of dimkey => Dimension to use in `dim` and `layerdims`
        dimdict = _dimdict(ds, crs, mappedcrs)
        refdims = isnokw(refdims) || isnothing(refdims) ? () : refdims
        metadata = isnokw(metadata) ? _metadata(ds) : metadata
        layerdims_vec = isnokw(layerdims) ? _layerdims(ds; layers, dimdict) : layerdims
        dims = _sort_by_layerdims(isnokw(dims) ? _dims(ds, dimdict) : dims, layerdims_vec)
        layermetadata_vec = if isnokw(layermetadata)
            _layermetadata(ds; layers)
        else
            layermetadata isa NamedTuple ? collect(layermetadata) : map(_ -> NoKW(), fn)
        end
        name = Tuple(map(Symbol, layers.names))
        NT = NamedTuple{name}
        layer_mvs = map(Rasters.missingval, layers.vars, layermetadata_vec)
        missingval_vec = if isnokw(missingval)
            layer_mvs .=> missing
        elseif missingval isa NamedTuple
            keys(missingval) == name || throw(ArgumentError("`missingval` names $(keys(missingval)) do not match layer names $name")) 
            layer_mvs .=> collect(missingval)
        elseif missingval === Rasters.missingval
            layer_mvs
        else
            layer_mvs .=> (missingval,) # Wrap in case its not iterable
        end::Vector
        eltype_vec = map(eltype, layers.vars)
        mod_vec = _stack_mods(eltype_vec, layermetadata_vec, missingval_vec; scaled, coerce)
        data = if lazy
            vars = ntuple(i -> layers.vars[i], length(name))
            mods = ntuple(i -> mod_vec[i], length(name))
            FileStack{typeof(source)}(ds, filename; name, group, mods, vars)
        else
            map(layers.vars, layermetadata_vec, mod_vec) do var, md, mod
                modvar = _maybe_modify(var, mod)
                checkmem && _checkobjmem(modvar)
                x = Array(modvar)
                x isa AbstractArray ? x : fill(x) # Catch an NCDatasets bug
            end |> NT
        end
        mv_outer = NT(map(_outer_missingval, mod_vec))
        return data, (; 
            dims,
            refdims,
            layerdims=NT(layerdims_vec),
            metadata,
            layermetadata=NT(layermetadata_vec),
            missingval=mv_outer,
        )
    end
    return RasterStack(data; field_kw..., kw...)
end


# Try to sort the dimensions by layer dimension into a sensible
# order that applies without permutation, preferencing the layers
# with most dimensions, and those that come first.
# Intentionally not type-stable
function _sort_by_layerdims(dims, layerdims)
    dimlist = union(layerdims)
    currentorder = nothing
    for i in length(dims):-1:1
        for ldims in dimlist
            length(ldims) == i || continue
            currentorder = _merge_dimorder(ldims, currentorder)
        end
    end
    return DD.dims(dims, currentorder)
end

_merge_dimorder(neworder, ::Nothing) = neworder
function _merge_dimorder(neworder, currentorder)
    ods = otherdims(neworder, currentorder)
    outorder = currentorder
    for od in ods
        # Get the dims position in current order
        i = findfirst(d -> d == od, neworder)
        found = false
        # Find the next dimension that is in the outorder
        for j in 1:length(ods)
            if length(neworder) >= (i + j)
                nextd = neworder[i + j]
                if nextd in outorder
                    n = dimnum(outorder, nextd)
                    outorder = (outorder[1:n-1]..., od, outorder[n:end]...)
                    found = true
                    break
                end
            end
        end
        if !found
            outorder = (outorder..., od)
        end
    end
    return outorder
end


function _layerkeysfromdim(A, dim)
    hasdim(A, dim) || throw(ArgumentError("`layersfrom` dim `$(name(dim))` not found in `$(map(basetypeof, dims(A)))`"))
    vals = parent(lookup(A, dim))
    l = length(vals)
    if l > MAX_STACK_SIZE
        D = basetypeof(dim)
        options = map(basetypeof, dims(otherdims(A, dim), d -> length(d) <= MAX_STACK_SIZE))
        toolong = "Lookup of `layersfrom=$D` is too long to use for stack layers: $l"
        if length(options) > 0
            throw(ArgumentError("$toolong. Choose a different dimension from: $options"))
        else
            throw(ArgumentError("$toolong. Maybe try a simple `Raster` ?"))
        end
    end
    map(vals) do x
        if x isa Number
            Symbol(string(DD.name(dim), "_", x))
        else
            Symbol(x)
        end
    end
end

Base.convert(::Type{RasterStack}, src::AbstractDimStack) = RasterStack(src)

Raster(stack::RasterStack) = cat(values(stack)...; dims=Band([keys(stack)...]))

defaultcrs(::Source, crs) = crs
defaultcrs(s::Source, ::NoKW) = defaultcrs(s)
defaultcrs(x) = nothing
defaultmappedcrs(::Source, crs) = crs
defaultmappedcrs(s::Source, ::NoKW) = defaultmappedcrs(s)
defaultmappedcrs(::Source) = nothing
