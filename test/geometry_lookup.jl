using NaturalEarth
using Test

using Rasters, DimensionalData
import Proj

@testset "construction" begin
    # fetch land polygons from Natural Earth
    land_polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry
    # create a DimVector from the land polygons
    gl = GeometryLookup(land_polygons)
    @test crs(gl) == EPSG(4326)
    @test all(GO.equals, zip(val(gl), land_polygons))
end

@testset "reprojecting a GeometryLookup" begin
    # fetch land polygons from Natural Earth
    land_polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry
    # create a DimVector from the land polygons
    dv = Raster(rand(Dim{:Geometry}(GeometryLookup(land_polygons))))
    # reproject the full vector data cube (vector data vector, in this case :D)
    target_crs = ProjString("+proj=wintri +type=crs")
    reprojected_via_rasters = val(dims(reproject(target_crs, dv), Dim{:Geometry}))
    reprojected_via_geometryops = GO.reproject(land_polygons; source_crs = EPSG(4326), target_crs = target_crs) 
    # test that the reprojected geometries are the same
    @test all(splat(GO.equals), zip(
        reprojected_via_rasters,  # reproject the vdc, get the geometries from it
        reprojected_via_geometryops # reproject the geometries directly
        )
    )
end
