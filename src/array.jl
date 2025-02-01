const CHECKMEM = Ref(true)

"""
    checkmem!(x::Bool)

Set `checkmem` to `true` or `false`.

In some architectures memory reporting may be wrong and you may
wish to disable memory checks.

This setting can be overridden with the `checkmem` keyword, where applicable.
"""
function checkmem!(checkmem::Bool)
    !checkmem || @warn "Setting `checkmem` to `false` globally may lead to out-of-memory errors or system crashes"
    CHECKMEM[] = checkmem
    return checkmem
end

const FLATTEN_SELECT = FileArray
const FLATTEN_IGNORE = Union{Dict,Set,Base.MultiplicativeInverses.SignedMultiplicativeInverse,Array}

"""
    AbstractRaster <: DimensionalData.AbstractDimArray

Abstract supertype for objects that wrap an array (or location of an array)
and metadata about its contents. It may be memory or hold a `FileArray`, which
holds the filename, and is only opened when required.

`AbstractRaster`s inherit from [`AbstractDimArray`]($DDarraydocs)
from DimensionalData.jl. They can be indexed as regular Julia arrays or with
DimensionalData.jl [`Dimension`]($DDdimdocs)s. They will plot as a heatmap in
Plots.jl with correct coordinates and labels, even after slicing with
`getindex` or `view`. `getindex` on a `AbstractRaster` will always return
a memory-backed `Raster`.
"""
abstract type AbstractRaster{T,N,D,A} <: AbstractDimArray{T,N,D,A} end

# Interface methods ###########################################################
"""
    missingval(x)

Returns the value representing missing data in the dataset
"""
function missingval end
missingval(_) = nothing
missingval(::AbstractArray{T}) where T = Missing <: T ? missing : nothing
missingval(A::AbstractRaster) = A.missingval

# The filename might be buried somewhere in a DiskArray wrapper, so try to get it
function filename(A::AbstractRaster)
    arrays = Flatten.flatten(parent(A), FLATTEN_SELECT, FLATTEN_IGNORE)
    if length(arrays) == 0
        nothing
    else
        # TODO its not clear which array supplies the filename if there are multiple in a broadcast
        filename(first(arrays))
    end
end
filename(::AbstractArray) = nothing # Fallback
filename(A::DiskArrays.AbstractDiskArray) = filename(parent(A))

cleanreturn(A::AbstractRaster) = rebuild(A, cleanreturn(parent(A)))
cleanreturn(x) = x

isdisk(A::AbstractRaster) = parent(A) isa DiskArrays.AbstractDiskArray
isdisk(x) = false
ismem(A::AbstractRaster) = !isdisk(A)

function Base.:(==)(A::AbstractRaster{T,N}, B::AbstractRaster{T,N}) where {T,N}
    size(A) == size(B) && all(A .== B)
end
for f in (:mappedbounds, :projectedbounds, :mappedindex, :projectedindex)
    @eval ($f)(A::AbstractRaster, dims_) = ($f)(dims(A, dims_))
    @eval ($f)(A::AbstractRaster) = ($f)(dims(A))
end

# DimensionalData methods

DD.units(A::AbstractRaster) = get(metadata(A), :units, nothing)
# Rebuild all types of AbstractRaster as Raster
function DD.rebuild(
    A::AbstractRaster, data, dims::Tuple, refdims, name,
    metadata, missingval=missingval(A)
)
    missingval1 = _fix_missingval(eltype(data), missingval)
    Raster(data, dims, refdims, name, metadata, missingval1)
end
function DD.rebuild(A::AbstractRaster;
    data=parent(A), dims=dims(A), refdims=refdims(A), name=name(A),
    metadata=metadata(A), missingval=missingval(A)
)
    rebuild(A, data, dims, refdims, name, metadata, missingval)
end

function DD.modify(f, A::AbstractRaster)
    newdata = if isdisk(A) # TODO may have to avoid calling `open` on DiskArray
        open(A) do O
            f(parent(O))
        end
    else
        f(parent(A))
    end
    size(newdata) == size(A) || error("$f returns an array with size $(size(newdata)) when the original size was $(size(A))")
    return rebuild(A, newdata)
end

DD.DimTable(As::Tuple{<:AbstractRaster,Vararg{AbstractRaster}}) =
    DD.DimTable(DimStack(map(read, As)))

# DiskArrays methods

DiskArrays.eachchunk(A::AbstractRaster) = DiskArrays.eachchunk(parent(A))
DiskArrays.haschunks(A::AbstractRaster) = DiskArrays.haschunks(parent(A))
DA.readblock!(A::AbstractRaster, dst, r::AbstractUnitRange...) =
    DA.readblock!(parent(A), dst, r...)
DA.writeblock!(A::AbstractRaster, src, r::AbstractUnitRange...) =
    DA.writeblock!(parent(A), src, r...)

# Base methods

Base.parent(A::AbstractRaster) = A.data
# Make sure a disk-based Raster is open before Array/collect
Base.Array(A::AbstractRaster) = open(O -> Array(parent(O)), A)
Base.collect(A::AbstractRaster) = open(O -> collect(parent(O)), A)

"""
    open(f, A::AbstractRaster; write=false)

`open` is used to open any `lazy=true` `AbstractRaster` and do multiple operations
on it in a safe way. The `write` keyword opens the file in write lookup so that it
can be altered on disk using e.g. a broadcast.

`f` is a method that accepts a single argument - an `Raster` object
which is just an `AbstractRaster` that holds an open disk-based object.
Often it will be a `do` block:

`lazy=false` (in-memory) rasters will ignore `open` and pass themselves to `f`.

```julia
# A is an `Raster` wrapping the opened disk-based object.
open(Raster(filepath); write=true) do A
    mask!(A; with=maskfile)
    A[I...] .*= 2
    # ...  other things you need to do with the open file
end
```

By using a do block to open files we ensure they are always closed again
after we finish working with them.
"""
function Base.open(f::Function, A::AbstractRaster; kw...)
    # Open FileArray to expose the actual dataset object, even inside nested wrappers
    fas = Flatten.flatten(parent(A), FLATTEN_SELECT, FLATTEN_IGNORE)
    if fas == ()
        f(Raster(parent(A), dims(A), refdims(A), name(A), metadata(A), missingval(A)))
    else
        if length(fas) == 1
            _open_one(f, A, fas[1]; kw...)
        else
            _open_many(f, A, fas; kw...)
        end
    end
end

function _open_one(f, A::AbstractRaster, fa::FileArray; kw...)
    open(fa; kw...) do x
        # Rewrap the opened object where the FileArray was nested in the parent array
        data = Flatten.reconstruct(parent(A), (x,), FLATTEN_SELECT, FLATTEN_IGNORE)
        openraster = Raster(data, dims(A), refdims(A), name(A), metadata(A), missingval(A))
        f(openraster)
    end
end

_open_many(f, A::AbstractRaster, fas::Tuple; kw...) = _open_many(f, A, fas, (); kw...)
function _open_many(f, A::AbstractRaster, fas::Tuple, oas::Tuple; kw...)
    open(fas[1]; kw...) do oa
        _open_many(f, A, Base.tail(fas), (oas..., oa); kw...)
    end
end
function _open_many(f, A::AbstractRaster, fas::Tuple{}, oas::Tuple; kw...)
    data = Flatten.reconstruct(parent(A), oas, FLATTEN_SELECT, FLATTEN_IGNORE)
    openraster = Raster(data, dims(A), refdims(A), name(A), metadata(A), missingval(A))
    f(openraster)
end

# Concrete implementation ######################################################

"""
    Raster <: AbstractRaster

    Raster(filepath::String; kw...)
    Raster(A::AbstractDimArray; kw...)
    Raster(A::AbstractArray, dims; kw...)

A generic [`AbstractRaster`](@ref) for spatial/raster array data. It can hold
either memory-backed arrays or, if `lazy=true`, a [`FileArray`](@ref),
which stores the `String` path to an unopened file.

If `lazy=true`, the file will only be opened lazily when it is indexed with `getindex`
or when `read(A)` is called. Broadcasting, taking a view, reversing, and most other
methods will _not_ load data from disk; they will be applied later, lazily.

# Arguments

- `dims`: `Tuple` of `Dimension`s needed when an `AbstractArray` is used.

# Keywords

$NAME_KEYWORD
$GROUP_KEYWORD
$MISSINGVAL_KEYWORD
$METADATA_KEYWORD
$CONSTRUCTOR_CRS_KEYWORD
$CONSTRUCTOR_MAPPEDCRS_KEYWORD
$REFDIMS_KEYWORD

When a filepath `String` is used:
$DROPBAND_KEYWORD
$LAZY_KEYWORD
$SOURCE_KEYWORD
$SCALED_KEYWORD
$RAW_KEYWORD

When A is an `AbstractDimArray`:
- `data`: can replace the data in an existing `AbstractRaster`
"""
struct Raster{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me,Mi<:Union{T,Nothing}} <: AbstractRaster{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
    function Raster(
        data::A, dims::D, refdims::R, name::Na, metadata::Me, missingval::Mi
    ) where {D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me,Mi} where {T,N}
        DD.checkdims(data, dims)
        missingval1 = _fix_missingval(T, missingval)
        new{T,N,D,R,A,Na,Me,typeof(missingval1)}(data, dims, refdims, name, metadata, missingval1)
    end
end
# Create a Raster from and AbstractArray and dims
function Raster(A::AbstractArray{T,N}, dims::Tuple;
    refdims=(),
    name=Symbol(""),
    metadata=NoMetadata(),
    missingval=Missing <: T ? missing : nothing,
    crs=nokw,
    mappedcrs=nokw
)::Raster{T,N} where {T,N}
    A = Raster(A, Dimensions.format(dims, A), refdims, name, metadata, missingval)
    A = isnokw(crs) ? A : setcrs(A, crs)
    A = isnokw(mappedcrs) ? A : setmappedcrs(A, mappedcrs)
    return A
end
# Create a Raster from and AbstractVector and dims,
# reshaping the Vector to match the dimensions
function Raster(A::AbstractArray{T,1}, dims::Tuple{<:Dimension,<:Dimension,Vararg};
    kw...
)::Raster{T,length(dims)} where T
    Raster(reshape(A, map(length, dims)), dims; kw...)
end
Raster(A::AbstractArray{<:Any,1}, dim::Dimension; kw...) = Raster(A, (dim,); kw...)
# Load a Raster from a table
function Raster(table, dims::Tuple;
    name=nokw,
    kw...
)
    Tables.istable(table) || throw(ArgumentError("First argument to `Raster` is not a table or other known object: $table"))
    name = isnokw(name) ? first(_not_a_dimcol(table, dims)) : name
    cols = Tables.columns(table)
    A = reshape(cols[name], map(length, dims))
    return Raster(A, dims; name, kw...)
end
# Load a Raster from another AbstractArray with `dims` as keyword
Raster(A::AbstractArray; dims, kw...) = Raster(A, dims; kw...)::Raster
# Load a Raster from another AbstractDimArray
function Raster(A::AbstractDimArray;
    data=parent(A),
    dims=dims(A),
    refdims=refdims(A),
    name=name(A),
    metadata=metadata(A),
    missingval=missingval(A),
    kw...
)::Raster
    return Raster(data, dims; refdims, name, metadata, missingval, kw...)
end
# Load a Raster from a string filename and predefined dimensions
function Raster(filename::AbstractString, dims::Tuple{<:Dimension,<:Dimension,Vararg};
    kw...
)::Raster
    Raster(filename; dims, kw...)
end
# Load a Raster from a string filename
function Raster(filename::AbstractString;
    source=nokw,
    kw...
)
    source = sourcetrait(filename, source)
    _open(filename; source, mod=NoMod()) do ds
        Raster(ds, filename; source, kw...)
    end::Raster
end
# Load a Raster from an opened Dataset
function Raster(ds, filename::AbstractString;
    dims=nokw,
    refdims=(),
    name=nokw,
    group=nokw,
    metadata=nokw,
    missingval=nokw,
    crs=nokw,
    mappedcrs=nokw,
    source=nokw,
    replace_missing=nokw,
    coerce=convert,
    scaled::Union{Bool,NoKW}=nokw,
    verbose::Bool=true,
    write::Bool=false,
    lazy::Bool=false,
    dropband::Bool=true,
    checkmem::Bool=CHECKMEM[],
    raw::Bool=false,
    mod=nokw,
    f=identity,
)::Raster
    _maybe_warn_replace_missing(replace_missing)
    # `raw` option will ignore `scaled` and `missingval`
    scaled, missingval = _raw_check(raw, scaled, missingval, verbose)
    # TODO use a clearer name for this
    name1 = filekey(ds, name)
    # Detect the source from filename
    source = sourcetrait(filename, source)
    # Open the dataset and variable specified by `name`, at `group` level if provided
    # At this level we do not apply `mod`.
    data_out, dims_out, metadata_out, missingval_out = _open(source, ds; name=name1, group, mod=NoMod()) do var
        metadata_out = isnokw(metadata) ? _metadata(var) : metadata
        missingval_out = _read_missingval_pair(var, metadata_out, missingval)
        # Generate mod for scaling
        mod = isnokw(mod) ? _mod(eltype(var), metadata_out, missingval_out; scaled, coerce) : mod
        # Define or load the data array
        data_out = if lazy
            # Define a lay FileArray
            FileArray{typeof(source)}(var, filename;
                name=name1, group, mod, write
            )
        else
            modvar = _maybe_modify(var, mod)
            # Check the data will fit in memory
            checkmem && _checkobjmem(modvar)
            # Move the modified array to memory
            Array(modvar)
        end
        # Generate dims
        dims_out = isnokw(dims) ? _dims(var, crs, mappedcrs) : format(dims, data_out)
        # Return the data to the parent function
        mv_outer = _outer_missingval(mod)
        data_out, dims_out, metadata_out, mv_outer
    end
    # Use name or an empty Symbol
    name_out = name1 isa Union{NoKW,Nothing} ? Symbol("") : Symbol(name1)
    # Define the raster
    raster = Raster(data_out, dims_out, refdims, name_out, metadata_out, missingval_out)
    # Maybe drop a single band dimension
    return _maybe_drop_single_band(raster, dropband, lazy)
end

filekey(ds, name) = name
filekey(filename::String) = Symbol(splitext(basename(filename))[1])

# Add a `dimconstructor` method so `AbstractProjected` lookups create a Raster
# TODO this should be unwrapped to `DD.lookupconstructor` to avoid future ambiguities
DD.dimconstructor(::Tuple{<:Dimension{<:AbstractProjected},Vararg{Dimension}}) = Raster