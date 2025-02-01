RA.sourcetrait(::ZD.ZarrVariable) = Zarrsource()
RA.sourcetrait(::ZD.ZarrDataset) = Zarrsource()
RA.sourceconstructor(::Type{Zarrsource}) = ZD.ZarrDataset
RA.checkfilename(::Zarrsource, filename) =
    isfile(filename) || isdir(filename) || RA._isurl(filename) || RA._filenotfound_error(filename)
# In ZarrDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{Zarrsource}) = nothing

# TODO move this upstream to CommonDataModel.jl and Datasets packages
function RA.missingval(md::RA.Metadata{<:Zarrsource})
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
        if !isnothing(mv) 
            fv == mv || @warn "Both '_FillValue' $fv and 'missing_value' $mv were found. Currently we only use the first."
        end
        return fv
    end
end