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
missingval(x) = missing
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
filename(A::AbstractArray) = nothing # Fallback

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
    Raster(data, dims, refdims, name, metadata, missingval)
end
function DD.rebuild(A::AbstractRaster;
    data=parent(A), dims=dims(A), refdims=refdims(A), name=name(A),
    metadata=metadata(A), missingval=missingval(A)
)
    rebuild(A, data, dims, refdims, name, metadata, missingval)
end

function DD.modify(f, A::AbstractRaster)
    newdata = if isdisk(A) 
        open(A) do O
            f(parent(O))
        end
    else
        f(parent(A))
    end
    size(newdata) == size(A) || error("$f returns an array with size $(size(newdata)) when the original size was $(size(A))")
    return rebuild(A, newdata)
end

function DD.DimTable(As::Tuple{<:AbstractRaster,Vararg{<:AbstractRaster}}...)
    DD.DimTable(DimStack(map(read, As...)))
end

# DiskArrays methods

DiskArrays.eachchunk(A::AbstractRaster) = DiskArrays.eachchunk(parent(A))
DiskArrays.haschunks(A::AbstractRaster) = DiskArrays.haschunks(parent(A))
function DA.readblock!(A::AbstractRaster, dst, r::AbstractUnitRange...)
    DA.readblock!(parent(A), dst, r...)
end
function DA.writeblock!(A::AbstractRaster, src, r::AbstractUnitRange...) 
    DA.writeblock!(parent(A), src, r...)
end

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
    Raster <: AbsractRaster

    Raster(filepath::AbstractString, dims; kw...)
    Raster(A::AbstractArray{T,N}, dims; kw...)
    Raster(A::AbstractRaster; kw...)

A generic [`AbstractRaster`](@ref) for spatial/raster array data. It may hold
memory-backed arrays or [`FileArray`](@ref), that simply holds the `String` path
to an unopened file. This will only be opened lazily when it is indexed with `getindex`
or when `read(A)` is called. Broadcasting, taking a view, reversing and most other 
methods _do not_ load data from disk: they are applied later, lazily.

# Keywords

- `dims`: `Tuple` of `Dimension`s for the array.
- `lazy`: A `Bool` specifying if to load the stack lazily from disk. `false` by default.
- `name`: `Symbol` name for the array, which will also retreive named layers if `Raster`
    is used on a multi-layered file like a NetCDF.
- `missingval`: value reprsenting missing data, normally detected form the file. Set manually
    when you know the value is not specified or is incorrect. This will *not* change any
    values in the raster, it simply assigns which value is treated as missing. To replace all of
    the missing values in the raster, use [`replace_missing`](@ref).
- `metadata`: `ArrayMetadata` object for the array, or `NoMetadata()`.
- `crs`: the coordinate reference system of  the objects `XDim`/`YDim` dimensions. 
    Only set this if you know the detected crs is incrorrect, or it is not present in
    the file. The `crs` is expected to be a GeoFormatTypes.jl `CRS` or `Mixed` `GeoFormat` type. 
- `mappedcrs`: the mapped coordinate reference system of the objects `XDim`/`YDim` dimensions.
    for `Mapped` lookups these are the actual values of the index. For `Projected` lookups
    this can be used to index in eg. `EPSG(4326)` lat/lon values, having it converted automatically.
    Only set this if the detected `mappedcrs` in incorrect, or the file does not have a `mappedcrs`,
    e.g. a tiff. The `mappedcrs` is expected to be a GeoFormatTypes.jl `CRS` or `Mixed` `GeoFormat` type. 
- `dropband`: drop single band dimensions. `true` by default.

# Internal Keywords

In some cases it is possible to set these keywords as well.

- `data`: can replace the data in an `AbstractRaster`
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from, defaulting to `()`.
"""
struct Raster{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me,Mi} <: AbstractRaster{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
end
function Raster(A::AbstractArray, dims::Tuple;
    refdims=(), name=Symbol(""), metadata=NoMetadata(), missingval=missing,
    crs=nothing, mappedcrs=nothing
)
    A = Raster(A, Dimensions.format(dims, A), refdims, name, metadata, missingval)
    A = isnothing(crs) ? A : setcrs(A, crs)
    A = isnothing(mappedcrs) ? A : setmappedcrs(A, mappedcrs)
    return A
end
function Raster(A::AbstractArray{<:Any,1}, dims::Tuple{<:Dimension,<:Dimension,Vararg}; kw...)
    Raster(reshape(A, map(length, dims)), dims; kw...)
end
function Raster(table, dims::Tuple; name=first(_not_a_dimcol(table, dims)), kw...)
    Tables.istable(table) || throw(ArgumentError("First argument to `Raster` is not a table or other known object: $table"))
    isnothing(name) && throw(UndefKeywordError(:name))
    cols = Tables.columns(table)
    A = reshape(cols[name], map(length, dims))
    return Raster(A, dims; name, kw...)
end
Raster(A::AbstractArray; dims, kw...) = Raster(A, dims; kw...)
function Raster(A::AbstractDimArray;
    data=parent(A), dims=dims(A), refdims=refdims(A),
    name=name(A), metadata=metadata(A), missingval=missingval(A), kw...
)
    return Raster(data, dims; refdims, name, metadata, missingval, kw...)
end
function Raster(filename::AbstractString, dims::Tuple{<:Dimension,<:Dimension,Vararg}; kw...)
    Raster(filename; dims, kw...)
end
function Raster(filename::AbstractString;
    name=nothing, key=name, source=nothing, kw...
)
    source = isnothing(source) ? _sourcetype(filename) : _sourcetype(source)
    _open(filename; source) do ds
        key = filekey(ds, key)
        Raster(ds, filename, key; kw...)
    end
end
function Raster(ds, filename::AbstractString, key=nothing;
    crs=nothing, mappedcrs=nothing, dims=nothing, refdims=(),
    name=Symbol(key isa Nothing ? "" : string(key)),
    metadata=metadata(ds), missingval=missingval(ds), 
    source=nothing, write=false, lazy=false, dropband=true,
)
    source = isnothing(source) ? _sourcetype(filename) : _sourcetype(source)
    crs = defaultcrs(source, crs)
    mappedcrs = defaultmappedcrs(source, mappedcrs)
    dims = dims isa Nothing ? DD.dims(ds, crs, mappedcrs) : dims
    data = if lazy 
        FileArray(ds, filename; key, write)
    else
        _open(Array, source, ds; key)
    end
    raster =  Raster(data, dims, refdims, name, metadata, missingval)
    return dropband ? _drop_single_band(raster, lazy) : raster
end

function Raster{T, N, D, R, A, Na, Me, Mi}(ras::Raster{T, N, D, R, A, Na, Me, Mi}) where {T, N, D, R, A, Na, Me, Mi}
    return Raster(ras.data, ras.dims, ras.refdims, ras.name, ras.metadata, ras.missingval)
end

filekey(ds, key) = key
filekey(filename::String) = Symbol(splitext(basename(filename))[1])

DD.dimconstructor(::Tuple{<:Dimension{<:AbstractProjected},Vararg{<:Dimension}}) = Raster


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
