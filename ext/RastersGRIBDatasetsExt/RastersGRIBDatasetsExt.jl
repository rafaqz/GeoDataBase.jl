module RastersGRIBDatasetsExt

using Rasters
using GRIBDatasets
using DiskArrays

using Rasters: GRIBsource

const RA = Rasters
const DA = DiskArrays
const GDS = GRIBDatasets

include("gribdatasets_source.jl")

end
