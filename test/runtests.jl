using Rasters, Test, Aqua, SafeTestsets

@testset "Aqua" begin
    Aqua.test_ambiguities(Rasters)
    Aqua.test_unbound_args(Rasters)
    Aqua.test_stale_deps(Rasters)
    Aqua.test_undefined_exports(Rasters)
    Aqua.test_project_extras(Rasters)
    # Aqua.test_deps_compat(Rasters) # This breaks GDAL downstream tests
end

@time @safetestset "extensions" begin include("extensions.jl") end
@time @safetestset "array" begin include("array.jl") end
@time @safetestset "stack" begin include("stack.jl") end
@time @safetestset "series" begin include("series.jl") end
@time @safetestset "create" begin include("create.jl") end
@time @safetestset "utils" begin include("utils.jl") end
@time @safetestset "set" begin include("set.jl") end
@time @safetestset "adapt" begin include("adapt.jl") end
@time @safetestset "methods" begin include("methods.jl") end
@time @safetestset "aggregate" begin include("aggregate.jl") end
@time @safetestset "mosaic" begin include("mosaic.jl") end
@time @safetestset "rasterize" begin include("rasterize.jl") end
@time @safetestset "extract" begin include("extract.jl") end
@time @safetestset "reproject" begin include("reproject.jl") end
@time @safetestset "warp" begin include("warp.jl") end
@time @safetestset "cellarea" begin include("cellarea.jl") end

@time @safetestset "sources" begin include("sources/sources.jl") end
@time @safetestset "commondatamodel" begin include("sources/commondatamodel.jl") end
@time @safetestset "ncdatasets" begin include("sources/ncdatasets.jl") end
@time @safetestset "zarrdatasets" begin include("sources/zarr.jl") end
if !Sys.iswindows()
    # GRIBDatasets doesn't work on Windows for now
    @time @safetestset "gribdatasets" begin include("sources/gribdatasets.jl") end
    # GDAL Environment vars need to be set manually for windows, so skip for now
    @time @safetestset "gdal" begin include("sources/gdal.jl") end
    @time @safetestset "grd" begin include("sources/grd.jl") end
end
# Only test RasterDataSources locally for now, because CI downloads keep breaking
if !haskey(ENV, "CI")
    @time @safetestset "rasterdatasources" begin include("sources/rasterdatasources.jl") end
end
@time @safetestset "plot recipes" begin include("plotrecipes.jl") end
@time @safetestset "resample" begin include("resample.jl") end

using ArchGDAL
using CoordinateTransformations
using GRIBDatasets
using Makie
using NCDatasets
using Proj
using RasterDataSources
using StatsBase
using ZarrDatasets

@testset "extensions Aqua" begin
    for name in [
        :RastersArchGDALExt,
        :RastersCoordinateTransformationsExt,
        :RastersGRIBDatasetsExt,
        # :RastersMakieExt, TODO lots of ambiguities
        :RastersNCDatasetsExt,
        :RastersProjExt,
        :RastersRasterDataSourcesExt,
        :RastersStatsBaseExt,
        :RastersZarrDatasetsExt,
    ]
        ext = Base.get_extension(Rasters, name)
        Aqua.test_ambiguities(ext)
        Aqua.test_unbound_args(ext)
        Aqua.test_undefined_exports(ext)
    end
end