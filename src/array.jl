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
const FLATTEN_IGNORE = Union{Dict,Set,Base.MultiplicativeInverses.SignedMultiplicativeInverse}

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
missingval(_) = missing
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
ismem(x) = !isdisk(x)

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
    # Have to avoid calling `open` on CFDiskArray
    newdata = if isdisk(A) && !(parent(A) isa CFDiskArray)
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

- `name`: a `Symbol` name for the array, which will also retrieve the, alphabetically first, 
    named layer if `Raster` is used on a multi-layered file like a NetCDF. 
    If instead `RasterStack` is used to read the multi-layered file, by default, all variables 
    will be added to the stack.
$GROUP_KEYWORD 
- `missingval`: value reprsenting missing data, normally detected from the file. Set manually
    when you know the value is not specified or is incorrect. This will *not* change any
    values in the raster, it simply assigns which value is treated as missing. To replace all of
    the missing values in the raster, use [`replace_missing`](@ref).
- `metadata`: `Dict` or `Metadata` object for the array, or `NoMetadata()`.
$CONSTRUCTOR_CRS_KEYWORD 
$CONSTRUCTOR_MAPPEDCRS_KEYWORD 
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from, defaulting to `()`.
    Usually not needed.

When a filepath `String` is used:
$DROPBAND_KEYWORD
$LAZY_KEYWORD
$REPLACE_MISSING_KEYWORD
$SOURCE_KEYWORD
- `write`: defines the default `write` keyword value when calling `open` on the Raster. `false` by default.
    Only makes sense to use when `lazy=true`.

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
function Raster(A::AbstractArray{T,1}, dims::Tuple{<:Dimension,<:Dimension,Vararg};
    kw...
)::Raster{T,length(dims)} where T
    Raster(reshape(A, map(length, dims)), dims; kw...)
end
Raster(A::AbstractArray{<:Any,1}, dim::Dimension; kw...) = Raster(A, (dim,); kw...)
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
Raster(A::AbstractArray; dims, kw...) = Raster(A, dims; kw...)::Raster
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
function Raster(filename::AbstractString, dims::Tuple{<:Dimension,<:Dimension,Vararg};
    kw...
)::Raster
    Raster(filename; dims, kw...)
end
Raster(::UndefInitializer, ext::Extents.Extent; kw...) =
    Raster{Float64}(undef, ext::Extents.Extent; kw...)
function Raster{T}(::UndefInitializer, ext::Extents.Extent; 
    size=nothing, res=nothing, crs=nothing, mappedcrs=nothing, sampling=Points(), closed=false, kw...
) where T
    dims = _extent2dims(ext; size, res, crs, mappedcrs, sampling, closed)
    Raster{T}(undef, dims; kw...)
end
function Raster(filename::AbstractString; 
    source=nokw,
    kw...
)
    source = _sourcetrait(filename, source)
    _open(filename; source) do ds
        Raster(ds, filename; source, kw...)
    end::Raster
end
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
    replace_missing=false,
    write=false,
    lazy=false,
    dropband=true,
    checkmem=CHECKMEM[],
)::Raster
    name1 = filekey(ds, name)
    source = _sourcetrait(filename, source)
    data1, dims1, metadata1, missingval1 = _open(source, ds; name=name1, group) do var
        metadata1 = isnokw(metadata) ? _metadata(var) : metadata
        missingval1 = _fix_missingval(var, missingval)
        rm = replace_missing && !isnothing(missingval1)
        missingval2 = rm ? missing : missingval1
        data = if lazy
            A = FileArray{typeof(source)}(var, filename; name=name1, group, write)
            rm ? _replace_missing(A, missingval1) : A
        else
            checkmem && _checkobjmem(var)
            x = Array(rm ? _replace_missing(var, missingval1) : var)
            x isa AbstractArray ? x : fill(x) # Catch an NCDatasets bug
        end
        dims1 = isnokw(dims) ? _dims(var, crs, mappedcrs) : format(dims, data)
        data, dims1, metadata1, missingval2
    end
    name2 = name1 isa Union{NoKW,Nothing} ? Symbol("") : Symbol(name1)
    raster = Raster(data1, dims1, refdims, name2, metadata1, missingval1)
    return dropband ? _drop_single_band(raster, lazy) : raster
end

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

function _replace_missing(A::AbstractArray{T}, missingval) where T
    repmissing(x) = isequal(x, missingval) ? missing : x
    return repmissing.(A)
end

filekey(ds, name) = name
filekey(filename::String) = Symbol(splitext(basename(filename))[1])

DD.dimconstructor(::Tuple{<:Dimension{<:AbstractProjected},Vararg{Dimension}}) = Raster


function _drop_single_band(raster, lazy::Bool)
    if hasdim(raster, Band()) && size(raster, Band()) < 2
         if lazy
             return view(raster, Band(1)) # TODO fix dropdims in DiskArrays
         else
             return dropdims(raster; dims=Band())
         end
    else
         return raster
    end
end
