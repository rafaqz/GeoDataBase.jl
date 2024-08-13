const NCDAllowedType = Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64,Char,String}

function RA._check_allowed_type(::RA.NCDsource, eltyp)
    eltyp <: NCDAllowedType || throw(ArgumentError("""
    Element type $eltyp cannot be written to NetCDF. Convert it to one of $(Base.uniontypes(NCDAllowedType)),
    usually by broadcasting the desired type constructor over the `Raster`, e.g. `newrast = Float32.(rast)`"))
    """
    ))
end
function Base.write(filename::AbstractString, ::NCDsource, A::AbstractRaster;
    append=false,
    force=false,
    kw...
)
    mode = if append
        isfile(filename) ? "a" : "c"
    else
        RA.check_can_write(filename, force)
        "c"
    end
    mode  = !isfile(filename) || !append ? "c" : "a";
    ds = NCD.Dataset(filename, mode; attrib=RA._attribdict(metadata(A)))
    try
        RA._writevar!(ds, A; kw...)
    finally
        close(ds)
    end
    return filename
end
function Base.write(filename::AbstractString, ::NCDsource, s::AbstractRasterStack;
    append=false,
    force=false,
    missingval=nokw,
    maskingval=nokw,
    kw...
)
    mode = if append
        isfile(filename) ? "a" : "c"
    else
        RA.check_can_write(filename, force)
        "c"
    end
    ds = NCD.Dataset(filename, mode; attrib=RA._attribdict(metadata(s)))

    maskingval = RA._stack_nt(s, isnokw(maskingval) ? Rasters.missingval(s) : maskingval)
    missingval = RA._stack_missingvals(s, isnokw(missingval) ? maskingval : missingval)
    try
        map(keys(s)) do k
            _writevar!(ds, s[k]; 
                missingval=missingval[k], 
                maskingval=maskingval[k], 
                kw...
            )
        end
    finally
        close(ds)
    end
    return filename
end

Base.close(os::RA.OpenStack{NCDsource}) = NCD.close(RA.dataset(os))

function RA.OpenStack(fs::RA.FileStack{NCDsource,K}) where K
    RA.OpenStack{NCDsource,K}(NCD.Dataset(RA.filename(fs)))
end

function RA._open(f, ::NCDsource, filename::AbstractString; write=false, kw...)
    isfile(filename) || RA._isurl(filename) || RA._filenotfound_error(filename)
    mode = write ? "a" : "r"
    NCD.Dataset(filename, mode) do ds
        RA._open(f, NCDsource(), ds; kw...)
    end
end

RA._sourcetrait(::NCD.Dataset) = NCDsource()
RA._sourcetrait(::NCD.Variable) = NCDsource()

@inline function RA.get_scale(metadata::Metadata{NCDsource}, scaled::Bool)
    scale = scaled ? get(metadata, "scale_factor", nothing) : nothing
    offset = scaled ? get(metadata, "add_offset", nothing) : nothing
    return scale, offset
end

RA.missingval(var::NCD.Variable, args...) = _mv(CDM.attribs(var))
RA.missingval(var::NCD.Variable, md::Metadata{<:NCDsource}) = _mv(md)

# TODO: handle multiple missing values
function _mv(md)
    fv = get(md, "_FillValue", nothing)
    mv = get(md, "missing_value", nothing)
    if isnothing(fv)
        if mv isa Vector
            length(mv) > 1 && @warn "'missing_value' $mv has multiple values. Currently we only uses the first."
            return first(mv)
        else
            return mv
        end
    else
        if isnothing(mv) 
            fv == mv || @warn "Both '_FillValue' $fv and 'missing_value' $mv were found. Currently we only use the first."
        end
        return fv
    end
end

# Add a var array to a dataset before writing it.
function _writevar!(ds::AbstractDataset, A::AbstractRaster{T,N};
    verbose=true,
    missingval=nokw,
    maskingval=nokw,
    metadata=nokw,
    chunks=nokw,
    chunksizes=RA._chunks_to_tuple(A, dims(A), chunks),
    scale=nokw,
    offset=nokw,
    coerce=convert,
    eltype=Missings.nonmissingtype(T),
    write=true,
    name=DD.name(A),
    options=nokw,
    driver=nokw,
    kw...
) where {T,N}
    eltype <: NCDAllowedType || throw(ArgumentError("""
       Element type $eltype cannot be written to NetCDF. Convert it to one of $(Base.uniontypes(NCDAllowedType)),
       usually by broadcasting the desired type constructor over the `Raster`, e.g. `newrast = Float32.(rast)`"))
       """
    ))
    _def_dim_var!(ds, A)
    metadata = if isnokw(metadata) 
        DD.metadata(A)
    elseif isnothing(metadata)
        NoMetadata()
    else
        metadata
    end

    maskingval = isnokw(maskingval) ? RA.missingval(A) : maskingval
    missingval = isnokw(missingval) ? RA.missingval(A) : missingval
    missingval = if ismissing(missingval) 
        # See if there is a missing value in metadata
        mv = _mv(metadata)
        # But only use it if its the right type
        mv isa eltype ? mv : RA._writeable_missing(eltype; verbose=true)
    else
        missingval
    end

    attrib = RA._attribdict(metadata)
    # Scale and offset
    scale = if isnokw(scale) || isnothing(scale)
        delete!(attrib, "scale_factor")
        nothing
    else
        attrib["scale_factor"] = scale
    end
    offset = if isnokw(offset) || isnothing(offset)
        delete!(attrib, "add_offset")
        nothing
    else
        attrib["add_offset"] = offset
    end

    mod = RA._writer_mod(eltype; missingval, maskingval, scale, offset, coerce)

    if !isnothing(mod.missingval)
        attrib["_FillValue"] = missingval
    end

    key = if isnokw(name) || string(name) == ""
        UNNAMED_NCD_FILE_KEY
    else
        string(name)
    end

    dimnames = lowercase.(string.(map(RA.name, dims(A))))
    var = NCD.defVar(ds, key, eltype, dimnames; attrib=attrib, chunksizes, kw...)

    if write
        # Write with a DiskArays.jl broadcast
        RA._maybe_modify(var.var, mod) .= A
    end

    return nothing
end

_def_dim_var!(ds::AbstractDataset, A) = map(d -> _def_dim_var!(ds, d), dims(A))
function _def_dim_var!(ds::AbstractDataset, dim::Dimension)
    dimname = lowercase(string(DD.name(dim)))
    haskey(ds.dim, dimname) && return nothing
    NCD.defDim(ds, dimname, length(dim))
    lookup(dim) isa NoLookup && return nothing

    # Shift index before conversion to Mapped
    dim = RA._cdmshiftlocus(dim)
    if dim isa Y || dim isa X
        dim = convertlookup(Mapped, dim)
    end
    # Attributes
    attrib = RA._attribdict(metadata(dim))
    RA._cdm_set_axis_attrib!(attrib, dim)
    # Bounds variables
    if sampling(dim) isa Intervals
        bounds = Dimensions.dim2boundsmatrix(dim)
        boundskey = get(metadata(dim), :bounds, string(dimname, "_bnds"))
        push!(attrib, "bounds" => boundskey)
        NCD.defVar(ds, boundskey, bounds, ("bnds", dimname))
    end
    NCD.defVar(ds, dimname, Vector(index(dim)), (dimname,); attrib=attrib)
    return nothing
end

# precompilation

# const _NCDVar = NCDatasets.CFVariable{Union{Missing, Float32}, 3, NCDatasets.Variable{Float32, 3, NCDatasets.NCDataset}, NCDatasets.Attributes{NCDatasets.NCDataset{Nothing}}, NamedTuple{(:fillvalue, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor), Tuple{Float32, Nothing, Nothing, Nothing, Nothing, Nothing}}}

# function _precompile(::Type{NCDsource})
#     ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

#     precompile(Rasters.FileArray, (_NCDVar, String))
#     precompile(layerkeys, (NCDatasets.NCDataset{Nothing},))
#     precompile(dims, (_NCDVar,Symbol))
#     precompile(dims, (_NCDVar,Symbol,Nothing,Nothing))
#     precompile(dims, (_NCDVar,Symbol,Nothing,EPSG))
#     precompile(dims, (_NCDVar,Symbol,EPSG,EPSG))
#     precompile(_firstkey, (NCDatasets.NCDataset{Nothing},))
#     precompile(_cdmdim, (NCDatasets.NCDataset{Nothing}, Symbol, Nothing, Nothing))
#     precompile(_cdmdim, (NCDatasets.NCDataset{Nothing}, Symbol, Nothing, EPSG))
#     precompile(_cdmdim, (NCDatasets.NCDataset{Nothing}, Symbol, EPSG, EPSG))
#     precompile(Raster, (NCDatasets.NCDataset{Nothing}, String, Nothing))
#     precompile(Raster, (NCDatasets.NCDataset{Nothing}, String, Symbol))
#     precompile(Raster, (_NCDVar, String, Symbol))

#     precompile(Raster, (String,))
# end

# _precompile(NCDsource)
