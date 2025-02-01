using Rasters
using DimensionalData
using DimensionalData.Lookups
using DimensionalData.Dimensions
using Dates
using ZarrDatasets
using ZarrDatasets.Zarr
using Rasters: FileArray, FileStack, Zarrsource, crs, bounds, name, trim

path = "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v3.0.2/esdc-16d-2.5deg-46x72x1440-3.0.2.zarr"

@testset "Zarr Raster" begin

zraster = Raster(path; name="air_temperature_2m")
lazyarray = Raster(path; lazy=true, name="air_temperature_2m")
eagerarray = Raster(path; lazy=false, name="air_temperature_2m")
@test_throws ArgumentError Raster("notafile.zarr/")

@testset "lazyness" begin
    # Eager is the default
    @test parent(zraster) isa Array
    @test parent(lazyarray) isa FileArray
    @test parent(eagerarray) isa Array
end

@testset "read" begin
    @time A = read(lazyarray);
    @test A isa Raster
    @test parent(A) isa Array
    A2 = copy(A) .= 0
    @time read!(zraster, A2);
    A3 = copy(A) .= 0
    @time read!(zraster, A3)
    @test all(A .=== A2) 
    @test all(A .=== A3)
end

@testset "array properties" begin
    @test name.(dims(zraster)) == (:X, :Y, :Ti)
    @test length(dims(zraster, X)) == 144
    @test index(zraster,X) == collect(-178.75:2.5:178.75)
    # TODO the spatial bounds are strange, because the data is point data
    # We should find a dataset that has actual intervals

    @test bounds(zraster) == (
        (-178.75, 178.75),
        (-88.75, 88.75),
        (DateTime("1979-01-09T00:00:00"), DateTime("2021-12-27T00:00:00")),
    )
end

@testset "dimensions" begin
    @test ndims(zraster) == 3
    @test length.(dims(zraster)) == (144, 72, 989)
    @test dims(zraster) isa Tuple{<:X,<:Y,<:Ti}
    @test refdims(zraster) == ()
    @test val.(span(zraster)) == (2.5, 2.5, (nothing, nothing))
    @test typeof(lookup(zraster)) <: Tuple{<:Mapped,<:Mapped,<:Sampled}
end

@testset "other fields" begin
    @test ismissing(missingval(zraster))
    @test metadata(zraster)["original_name"] == "t2m"
    @test metadata(zraster) isa Rasters.Metadata{<:Rasters.CDMsource, Dict{String, Any}}
    @test name(zraster) == :air_temperature_2m
end

@testset "indexing" begin
    @test zraster[Ti(1)] isa Raster{<:Any,2}
    @test zraster[Y(1), Ti(1)] isa Raster{<:Any,1}
    @test zraster[X(1), Ti(1)] isa Raster{<:Any,1}
    @test zraster[X(1), Y(1), Ti(1)] == -28.866226f0 == parent(zraster)[1,1,1]
    @test zraster[X(30), Y(30), Ti(1)] isa Float32
    # Alaska
    @test zraster[Y(Near(-88.75)), X(Near(-178.74)), Ti(1)] ==-28.866226f0
    @test zraster[Ti(At(DateTime(1979,1,9))), X(At(-178.75)), Y(At(-88.75))] == -28.866226f0
end

@testset "CF conventions" begin
    path = tempname() * ".zarr"
    global_attrs = Dict(
        "Conventions" => "CF-1.7",
        "_FillValue" => -999.0,
        "units" => "K",
        "project" => "GOES",
    )
    zg = zgroup(path; attrs = global_attrs)
    x_attrs = Dict(
        "units"             => "rad",
        "add_offset"        => 0,
        "_ARRAY_DIMENSIONS" => ["x"],
        "axis"              => "X",
        "long_name"         => "GOES fixed grid projection x-coordinate",
        "scale_factor"      => 1,
        "_Netcdf4Dimid"     => 1,
        "standard_name"     => "projection_x_coordinate",
        "NAME"              => "x",
    )
    xs = zcreate(Float64, zg, "x", 100; attrs = x_attrs)
    xs .= LinRange(-5.434177815866809e6, 5.322929839864517e6, 100)

    y_attrs = Dict(
        "units"             => "rad",
        "add_offset"        => 5.258454335331338e6,
        "_ARRAY_DIMENSIONS" => ["y"],
        "axis"              => "Y",
        "long_name"         => "GOES fixed grid projection y-coordinate",
        "scale_factor"      => -1.0625617981243342e7,
        "_Netcdf4Dimid"     => 1,
        "standard_name"     => "projection_y_coordinate",
        "NAME"              => "y",
    )
    ys = zcreate(Float64, zg, "y", 100; attrs = y_attrs)
    ys .= LinRange(0, 1, 100)

    val_attrs = Dict(
        "_Netcdf4Dimid" => 0,
        "scale_factor" => 2.0,
        "add_offset" => 180.0,
        "_FillValue" => -1, 
        "units" => "K", 
        "_ARRAY_DIMENSIONS" => ["y", "x"], 
        "grid_mapping" => "goes_imager_projection",
        "valid_range" => [0, -6],
    )
    vals = zcreate(Float64, zg, "values", 100, 100; attrs = val_attrs)
    vals .= (data = rand(100, 100))
    vals[1, 1] = 1.0
    vals[end, end] = 0.0

    ra = Raster(path)

    @test extrema(ra) == (180.0, 182.0) # test scale and offset
    @test Rasters.GeoInterface.extent(ra) == Rasters.GeoInterface.Extent(X = (-5.434177815866809e6, 5.322929839864517e6), Y = (-5.367163645912004e6, 5.258454335331338e6))
    @test Rasters.isreverse(ra.dims[2])
    @test Rasters.isforward(ra.dims[1])
    @test extrema(ra.dims[1]) == extrema(xs)
    @test extrema(ra.dims[2]) == reverse(extrema(ys)) .* y_attrs["scale_factor"] .+ y_attrs["add_offset"]

    @testset "write" begin
        fn = tempname() * ".zarr"
        write(fn, ra)
        @test all(Raster(fn) .=== ra)
    end
end


@testset "RasterStack" begin
    st = RasterStack(path; lazy=true)
    @test all(st.snow_sublimation .=== Raster(path; name=:snow_sublimation))
end

end
