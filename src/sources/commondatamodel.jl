const CDM = CommonDataModel

const UNNAMED_FILE_KEY = "unnamed"

const CDM_DIM_MAP = Dict(
    "lat" => Y,
    "latitude" => Y,
    "lon" => X,
    "long" => X,
    "longitude" => X,
    "time" => Ti,
    "lev" => Z,
    "mlev" => Z,
    "level" => Z,
    "vertical" => Z,
    "x" => X,
    "y" => Y,
    "z" => Z,
    "band" => Band,
)

const CDM_AXIS_MAP = Dict(
    "X" => X,
    "Y" => Y,
    "Z" => Z,
    "T" => Ti,
)

const CDM_STANDARD_NAME_MAP = Dict(
    "longitude" => X,
    "latitude" => Y,
    "depth" => Z,
    "time" => Ti,
)

sourcetrait(var::CDM.CFVariable) = sourcetrait(var.var)
sourceconstructor(source::Source) = sourceconstructor(typeof(source))
function checkfilename end
function checkwritemode(::CDMsource, filename, append::Bool, force::Bool)
    if append
        isfile(filename) ? "a" : "c"
    else
        check_can_write(filename, force)
        "c"
    end
end

openmode(write::Bool) = write ? "a" : "r"

missingval(var::CDM.AbstractVariable, md::Metadata{<:CDMsource}) =
    missingval(md)
missingval(var::CDM.AbstractVariable, args...) = 
    missingval(Metadata{soucetrait(var)}(CDM.attribs(var)))

@inline function get_scale(metadata::Metadata{<:CDMsource}, scaled::Bool)
    scale = scaled ? get(metadata, "scale_factor", nothing) : nothing
    offset = scaled ? get(metadata, "add_offset", nothing) : nothing
    return scale, offset
end

# Rasters methods for CDM types ###############################

function FileStack{source}(ds::AbstractDataset, filename::AbstractString;
    write::Bool=false, 
    group=nokw,
    name::NTuple{N,Symbol}, 
    mods,
    vars,
) where {source<:CDMsource,N}
    T = NamedTuple{name,Tuple{map(_mod_eltype, vars, mods)...}}
    layersizes = map(size, vars)
    eachchunk = map(DiskArrays.eachchunk, vars)
    haschunks = map(DiskArrays.haschunks, vars)
    group = isnokw(group) ? nothing : group
    return FileStack{source,name,T}(filename, layersizes, group, eachchunk, haschunks, mods, write)
end

OpenStack(fs::FileStack{Source,K}) where {K,Source<:CDMsource} =
    OpenStack{Source,K}(sourceconstructor(Source())(filename(fs)), fs.mods)

function _open(f, source::CDMsource, filename::AbstractString; write=false, kw...)
    checkfilename(source, filename)
    ds = sourceconstructor(source)(filename)
    _open(f, source, ds; kw...)
end
function _open(f, source::CDMsource, ds::AbstractDataset; 
    name=nokw, 
    group=nothing, 
    mod=NoMod(), 
    kw...
)
    g = _getgroup(ds, group)
    x = if isnokw(name)
        cleanreturn(f(g)) 
    else
        v = CDM.variable(g, string(_name_or_firstname(g, name)))
        _open(f, source, v; mod)
    end
end
_open(f, ::CDMsource, var::AbstractArray; mod=NoMod(), kw...) = 
    cleanreturn(f(_maybe_modify(var, mod)))

# This allows arbitrary group nesting
_getgroup(ds, ::Union{Nothing,NoKW}) = ds
_getgroup(ds, group::Union{Symbol,AbstractString}) = ds.group[String(group)]
_getgroup(ds, group::Pair) = _getgroup(ds.group[String(group[1])], group[2])

filekey(ds::AbstractDataset, name::Union{String,Symbol}) = Symbol(name)
filekey(ds::AbstractDataset, name) = _name_or_firstname(ds, name)
cleanreturn(A::AbstractVariable) = Array(A)
haslayers(::CDMsource) = true
defaultcrs(::CDMsource) = EPSG(4326)
defaultmappedcrs(::CDMsource) = EPSG(4326)

function _nondimnames(ds)
    dimnames = CDM.dimnames(ds)
    toremove = if "bnds" in dimnames
        dimnames = setdiff(dimnames, ("bnds",))
        boundsnames = String[]
        for k in dimnames
            var = ds[k]
            attr = CDM.attribs(var)
            if haskey(attr, "bounds")
                push!(boundsnames, attr["bounds"])
            end
        end
        union(dimnames, boundsnames)::Vector{String}
    else
        collect(dimnames)::Vector{String}
    end
    # Maybe this should be fixed in ZarrDatasets but it works with this patch.
    nondim = collect(setdiff(keys(ds), toremove))
    return nondim
end


function _layers(ds::AbstractDataset, ::NoKW=nokw, ::NoKW=nokw)
    nondim = _nondimnames(ds)
    grid_mapping = String[]
    vars = map(k -> CDM.variable(ds, k), nondim)
    attrs = map(CDM.attribs, vars)
    for attr in attrs
        if haskey(attr, "grid_mapping")
            push!(grid_mapping, attr["grid_mapping"])
        end
    end
    bitinds = map(!in(grid_mapping), nondim)
    (;
        names=nondim[bitinds],
        vars=vars[bitinds],
        attrs=attrs[bitinds],
    )
end
_layers(ds::AbstractDataset, names, ::NoKW) = 
    _layers(ds, collect(names), nokw)
function _layers(ds::AbstractDataset, names::Vector, ::NoKW)
    vars = map(n -> CDM.variable(ds, n), names)
    attrs = map(CDM.attribs, vars)
    (; names, vars, attrs)
end
function _layers(ds::AbstractDataset, names, group)
    _layers(ds.group[group], names, nokw)
end

function _dims(var::AbstractVariable{<:Any,N}, crs=nokw, mappedcrs=nokw) where N
    dimnames = CDM.dimnames(var)
    ntuple(Val(N)) do i
        _cdmdim(CDM.dataset(var), dimnames[i], crs, mappedcrs)
    end
end
_metadata(var::AbstractVariable; attr=CDM.attribs(var)) =
    _metadatadict(sourcetrait(var), attr)

function _dimdict(ds::AbstractDataset, crs=nokw, mappedcrs=nokw)
    dimdict = Dict{String,Dimension}()
    for dimname in CDM.dimnames(ds)
        dimdict[dimname] = _cdmdim(ds, dimname, crs, mappedcrs)
    end
    return dimdict
end
function _dims(ds::AbstractDataset, dimdict::Dict)
    map(CDM.dimnames(ds)) do dimname
        dimdict[dimname]
    end |> Tuple
end
_metadata(ds::AbstractDataset; attr=CDM.attribs(ds)) =
    _metadatadict(sourcetrait(ds), attr)
function _layerdims(ds::AbstractDataset; layers, dimdict)
    map(layers.vars) do var
        map(CDM.dimnames(var)) do dimname
            basedims(dimdict[dimname])
        end |> Tuple
    end
end
function _layermetadata(ds::AbstractDataset; layers)
    map(layers.attrs) do attr
        md = _metadatadict(sourcetrait(ds), attr)
        if haskey(attr, "grid_mapping")
            if haskey(ds, attr["grid_mapping"])
                md["grid_mapping"] = Dict(CDM.attribs(ds[attr["grid_mapping"]]))
            else
                global_attrs = CDM.attribs(ds)
                if haskey(global_attrs, attr["grid_mapping"])
                    md["grid_mapping"] = global_attrs["grid_mapping"]
                end
            end     
        end
        md
    end
end


# Utils ########################################################################

# TODO don't load all keys here with _layers
_name_or_firstname(ds::AbstractDataset, name) = Symbol(name)
function _name_or_firstname(ds::AbstractDataset, name::Union{Nothing,NoKW}=nokw)
    names = _nondimnames(ds)
    if length(names) > 0
        return Symbol(first(names))
    else
        throw(ArgumentError("No non-dimension layers found in dataset with keys: $(keys(ds))"))
    end
end

function _cdmdim(ds, dimname::Key, crs=nokw, mappedcrs=nokw)
    if haskey(ds, dimname)
        var = ds[dimname]
        D = _cdmdimtype(CDM.attribs(var), dimname)
        lookup = _cdmlookup(ds, dimname, D, crs, mappedcrs)
        return D(lookup)
    else
        # The var doesn't exist. Maybe its `complex` or some other marker,
        # so make it a custom `Dim` with `NoLookup`
        len = _cdmfinddimlen(ds, dimname)
        len === nothing && _unuseddimerror()
        lookup = NoLookup(Base.OneTo(len))
        D = _cdmdimtype(NoMetadata(), dimname)
        return D(lookup)
    end
end

function _cdmfinddimlen(ds, dimname)
    for name in keys(ds)
        var = ds[name]
        dimnames = CDM.dimnames(var)
        if dimname in dimnames
            return size(var)[findfirst(==(dimname), dimnames)]
        end
    end
    return nothing
end

# Find the matching dimension constructor. If its an unknown name
# use the generic Dim with the dim name as type parameter
function _cdmdimtype(attrib, dimname)
    if haskey(attrib, "axis")
        k = attrib["axis"]
        if haskey(CDM_AXIS_MAP, k)
            return CDM_AXIS_MAP[k]
        end
    end
    if haskey(attrib, "standard_name")
        k = attrib["standard_name"]
        if haskey(CDM_STANDARD_NAME_MAP, k)
            return CDM_STANDARD_NAME_MAP[k]
        end
    end
    if haskey(CDM_DIM_MAP, dimname)
        return CDM_DIM_MAP[dimname]
    end
    return DD.basetypeof(DD.name2dim(Symbol(dimname)))
end

# _cdmlookup
# Generate a `Lookup` from a nCDM dim.
function _cdmlookup(ds::AbstractDataset, dimname, D::Type, crs, mappedcrs)
    var = ds[dimname]
    index = Missings.disallowmissing(var[:])
    attr = CDM.attribs(var)
    metadata = _metadatadict(sourcetrait(ds), attr)
    return _cdmlookup(ds, var, attr, dimname, D, index, metadata, crs, mappedcrs)
end
# For unknown types we just make a Categorical lookup
function _cdmlookup(ds::AbstractDataset, var, attr, dimname, D::Type, index::AbstractArray, metadata, crs, mappedcrs)
    Categorical(index; order=Unordered(), metadata=metadata)
end
# For Number and AbstractTime we generate order/span/sampling
# We need to include `Missing` in unions in case `_FillValue` is used
# on coordinate variables in a file and propagates here.
function _cdmlookup(
    ds::AbstractDataset, var, attr, dimname,
    D::Type, index::AbstractArray{<:Union{Missing,Number,Dates.AbstractTime}},
    metadata, crs, mappedcrs
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    order = LA.orderof(index)
    var = ds[dimname]
    # Detect lat/lon
    span, sampling = if eltype(index) <: Union{Missing,Dates.AbstractTime}
        _cdmperiod(index, metadata)
    else
        _cdmspan(index, order)
    end
    # We only use Explicit if the span is not Regular
    # This is important for things like rasterizatin and conversion 
    # to gdal to be easy, and selectors are faster.
    # TODO are there any possible floating point errors from this?
    if haskey(CDM.attribs(var), "bounds")
        span, sampling = if isregular(span)
            span, Intervals(Center())
        else
            boundskey = var.attrib["bounds"]
            boundsmatrix = Array(ds[boundskey])
            locus = if mapreduce(==, &, view(boundsmatrix, 1, :), index)
                Start()
            elseif mapreduce(==, &, view(boundsmatrix, 2, :), index)
                End()
            else
                Center()
            end
            Explicit(boundsmatrix), Intervals(locus)
        end
    end

    # We cant yet check CF standards crs, but we can at least check for units in lat/lon 
    if isnokw(mappedcrs) && get(metadata, "units", "") in ("degrees_north", "degrees_east")
        mappedcrs = EPSG(4326)
    end
    # Additionally, crs and mappedcrs should be identical for Regular lookups
    if isnokw(crs) && span isa Regular
        crs = mappedcrs
    end
    return _cdmlookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
end
# For X and Y use a Mapped <: AbstractSampled lookup
function _cdmlookup(
    D::Type{<:Union{<:XDim,<:YDim}}, index, order::Order, span, sampling, metadata, crs, mappedcrs
)
    # If the index is regularly spaced and there is no crs
    # then there is probably just one crs - the mappedcrs
    crs = if isnokw(crs) && span isa Regular
        mappedcrs
    else
        crs
    end
    dim = DD.basetypeof(D)()
    return Mapped(index; order, span, sampling, metadata, crs, mappedcrs, dim)
end
# Band dims have a Categorical lookup, with order
function _cdmlookup(D::Type{<:Band}, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Categorical(index, order, metadata)
end
# Otherwise use a regular Sampled lookup
function _cdmlookup(D::Type, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Sampled(index, order, span, sampling, metadata)
end

function _cdmspan(index, order)
    # Handle a length 1 index
    length(index) == 1 && return Regular(zero(eltype(index))), Points()

    step = if eltype(index) <: AbstractFloat
        # Calculate step, avoiding as many floating point errors as possible
        st = Base.step(Base.range(Float64(first(index)), Float64(last(index)); length = length(index)))
        st_rd = round(st, digits = Base.floor(Int,-log10(eps(eltype(index))))) # round to nearest digit within machine epsilon
        isapprox(st_rd, st; atol = eps(eltype(index))) ? st_rd : st # keep the rounded number if it is very close to the original
    else
        index[2] - index[1]
    end
    for i in 2:length(index)-1
        # If any step sizes don't match, its Irregular
        if !(index[i+1] - index[i] ≈ step)
            bounds = if length(index) > 1
                beginhalfcell = abs((index[2] - index[1]) * 0.5)
                endhalfcell = abs((index[end] - index[end-1]) * 0.5)
                if LA.isrev(order)
                    index[end] - endhalfcell, index[1] + beginhalfcell
                else
                    index[1] - beginhalfcell, index[end] + endhalfcell
                end
            else
                index[1], index[1]
            end
            return Irregular(bounds), Points()
        end
    end
    # Otherwise regular
    return Regular(step), Points()
end

# delta_t and ave_period are not CF standards, but CDC
function _cdmperiod(index, metadata::Metadata{<:CDMsource})
    if haskey(metadata, "delta_t")
        period = _parse_period(metadata["delta_t"])
        period isa Nothing || return Regular(period), Points()
    elseif haskey(metadata, "avg_period")
        period = _parse_period(metadata["avg_period"])
        period isa Nothing || return Regular(period), Intervals(Center())
    end
    return sampling = Irregular((nothing, nothing)), Points()
end

function _parse_period(period_str::String)
    regex = r"(\d\d\d\d)-(\d\d)-(\d\d) (\d\d):(\d\d):(\d\d)"
    mtch = match(regex, period_str)
    if mtch === nothing
        return nothing
    else
        vals = map(x -> parse(Int, x), mtch.captures)
        if length(vals) == 6
            y = Year(vals[1])
            mo = Month(vals[2])
            d = Day(vals[3])
            h = Hour(vals[4])
            mi = Minute(vals[5])
            s = Second(vals[6])
            compound = sum((y, mo, d, h, mi, s))
            if length(compound.periods) == 0
                return nothing
            elseif length(compound.periods) == 1
                return compound.periods[1]
            else
                return compound
            end
        else
            return nothing
        end
    end
end

function _attribdict(md::Metadata{<:CDMsource}) 
    attrib = Dict{String,Any}()
    for (k, v) in md
        v isa Tuple && continue
        attrib[string(k)] = v
    end
    return attrib
end
_attribdict(md) = Dict{String,Any}()

# Add axis and standard name attributes to dimension variables
# We need to get better at guaranteeing if X/Y is actually measured in `longitude/latitude`
# CF standards requires that we specify "units" if we use these standard names
_cdm_set_axis_attrib!(atr, dim::X) = atr["axis"] = "X" # at["standard_name"] = "longitude";
_cdm_set_axis_attrib!(atr, dim::Y) = atr["axis"] = "Y" # at["standard_name"] = "latitude";
_cdm_set_axis_attrib!(atr, dim::Z) = (atr["axis"] = "Z"; atr["standard_name"] = "depth")
_cdm_set_axis_attrib!(atr, dim::Ti) = (atr["axis"] = "T"; atr["standard_name"] = "time")
_cdm_set_axis_attrib!(atr, dim) = nothing

_cdmshiftlocus(dim::Dimension) = _cdmshiftlocus(lookup(dim), dim)
_cdmshiftlocus(::Lookup, dim::Dimension) = dim
function _cdmshiftlocus(lookup::AbstractSampled, dim::Dimension)
    # TODO improve this
    if span(lookup) isa Regular && sampling(lookup) isa Intervals
        # We cant easily shift a DateTime value
        if eltype(dim) isa Dates.AbstractDateTime
            if !(locus(dim) isa Center)
                @warn "To save to netcdf, DateTime values should be the interval Center, rather than the $(nameof(typeof(locus(dim))))"
            end
            dim
        else
            shiftlocus(Center(), dim)
        end
    else
        dim
    end
end

_unuseddimerror(dimname) = error("Dataset contains unused dimension $dimname")

function Base.write(filename::AbstractString, source::CDMsource, A::AbstractRaster;
    append=false,
    force=false,
    kw...
)
    mode = checkwritemode(source, filename, append, force)
    ds = sourceconstructor(source)(filename, mode; attrib=_attribdict(metadata(A)))
    try
        writevar!(ds, source, A; kw...)
    finally
        close(ds)
    end
    return filename
end
function Base.write(filename::AbstractString, source::Source, s::AbstractRasterStack{K,T};
    append=false,
    force=false,
    missingval=nokw,
    f=identity,
    kw...
) where {Source<:CDMsource,K,T}
    mode = checkwritemode(source, filename, append, force)
    ds = sourceconstructor(source)(filename, mode; attrib=_attribdict(metadata(s)))
    missingval = _stack_nt(s, isnokw(missingval) ? Rasters.missingval(s) : missingval)
    try
        map(keys(s)) do k
            writevar!(ds, source, s[k]; missingval=missingval[k], kw...)
        end
        f(OpenStack{Source,K,T}(ds))
    finally
        close(ds)
    end
    return filename
end

# Add a var array to a dataset before writing it.
function writevar!(ds::AbstractDataset, source::CDMsource, A::AbstractRaster{T,N};
    verbose=true,
    missingval=nokw,
    metadata=nokw,
    chunks=nokw,
    chunksizes=_chunks_to_tuple(A, dims(A), chunks),
    scale=nokw,
    offset=nokw,
    coerce=convert,
    eltype=Missings.nonmissingtype(T),
    write=true,
    name=DD.name(A),
    options=nokw,
    driver=nokw,
    f=identity,
    kw...
) where {T,N}
    _check_allowed_type(source, eltype)
    write = f === identity ? write : true
    _def_dim_var!(ds, A)
    metadata = if isnokw(metadata) 
        DD.metadata(A)
    elseif isnothing(metadata)
        NoMetadata()
    else
        metadata
    end

    missingval_pair = _write_missingval_pair(A, missingval; eltype, verbose, metadata)

    attrib = _attribdict(metadata)
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

    mod = _mod(eltype, missingval_pair, scale, offset, coerce)

    if !isnothing(missingval_pair[1])
        attrib["_FillValue"] = missingval_pair[1]
    end

    key = if isnokw(name) || string(name) == ""
        UNNAMED_FILE_KEY
    else
        string(name)
    end

    dimnames = lowercase.(string.(map(Rasters.name, dims(A))))
    var = CDM.defVar(ds, key, eltype, dimnames; attrib, chunksizes, kw...)

    if write
        m = _maybe_modify(var.var, mod)
        # Write with a DiskArays.jl broadcast
        m .= A
        # Apply `f` while the variable is open
        f(m)
    end

    return nothing
end

const CDMallowedType = Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64,Char,String}

_check_allowed_type(trait, eltyp) = nothing
function _check_allowed_type(::CDMsource, eltyp)
    eltyp <: CDMallowedType || throw(ArgumentError("""
    Element type $eltyp cannot be written to NetCDF. Convert it to one of $(Base.uniontypes(CDMallowedType)),
    usually by broadcasting the desired type constructor over the `Raster`, e.g. `newrast = Float32.(rast)`"))
    """
    ))
end


_def_dim_var!(ds::AbstractDataset, A) = map(d -> _def_dim_var!(ds, d), dims(A))
function _def_dim_var!(ds::AbstractDataset, dim::Dimension)
    dimname = lowercase(string(DD.name(dim)))
    haskey(ds.dim, dimname) && return nothing
    CDM.defDim(ds, dimname, length(dim))
    lookup(dim) isa NoLookup && return nothing

    # Shift index before conversion to Mapped
    dim = _cdmshiftlocus(dim)
    if dim isa Y || dim isa X
        dim = convertlookup(Mapped, dim)
    end
    # Attributes
    attrib = _attribdict(metadata(dim))
    _cdm_set_axis_attrib!(attrib, dim)
    # Bounds variables
    if sampling(dim) isa Intervals
        bounds = Dimensions.dim2boundsmatrix(dim)
        boundskey = get(metadata(dim), :bounds, string(dimname, "_bnds"))
        push!(attrib, "bounds" => boundskey)
        CDM.defVar(ds, boundskey, bounds, ("bnds", dimname))
    end
    CDM.defVar(ds, dimname, Vector(index(dim)), (dimname,); attrib=attrib)
    return nothing
end
