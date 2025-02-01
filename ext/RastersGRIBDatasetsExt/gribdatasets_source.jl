RA.sourcetrait(::GDS.Variable) = GRIBsource()
RA.sourcetrait(::GDS.GRIBDataset) = GRIBsource()

RA.sourceconstructor(::Type{GRIBsource}) = GDS.GRIBDataset

RA.checkwritemode(::GRIBsource, filename, append::Bool, force::Bool) =
    throw(ArgumentError("GRIBDatasets.jl does not support writing"))
RA.checkfilename(::GRIBsource, filename) =
    isfile(filename) || _filenotfound_error(filename)

# In GRIBDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{GRIBsource}) = nothing

RA.missingval(var::GDS.Variable, ::RA.Metadata{GRIBsource}) = _missingval(var)
RA.missingval(var::GDS.Variable, args...) = _missingval(var)

function _missingval(var::GDS.Variable{T}) where T
    mv = GDS.missing_value(var)
    T1 = promote_type(typeof(mv), T)
    return T1(mv)
end