using Documenter, GeoData, HDF5, ArchGDAL, NCDatasets

makedocs(
    modules = [GeoData],
    sitename = "GeoData.jl",
    strict = true,
)

deploydocs(
    repo = "github.com/rafaqz/GeoData.jl.git",
)
