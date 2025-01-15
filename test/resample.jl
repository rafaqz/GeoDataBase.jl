using Rasters, ArchGDAL, GeoInterface
using Extents
using Test
using Rasters.Lookups
import DimensionalData: @dim, YDim
@dim Lat YDim "Lat"

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

@testset "resample" begin
    raster_path = maybedownload("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")

    output_res = 0.0027
    output_crs = EPSG(4326)
    method = "near"

    ## Resample cea.tif manually with ArchGDAL
    wkt = convert(String, convert(WellKnownText, output_crs))
    AG_output = ArchGDAL.read(raster_path) do dataset
        ArchGDAL.gdalwarp([dataset], ["-t_srs", "$(wkt)",
                                "-tr", "$(output_res)", "$(output_res)",
                                "-r", "$(method)"]) do warped
            ArchGDAL.read(ArchGDAL.getband(warped, 1))
        end
    end

    # Resample cea.tif using resample
    cea = Raster(raster_path; missingval=0x00, name = :cea)
    raster_output = resample(cea; res=output_res, crs=output_crs, method)
    disk_output = resample(cea; res=output_res, crs=output_crs, method, filename="resample.tif")

    cea_permuted = permutedims(Raster(raster_path), (Y, X))
    permuted_output = resample(cea_permuted, output_res; crs=output_crs, method)

    # Compare ArchGDAL, resample and permuted resample 
    @test AG_output ==
        raster_output[Band(1)] ==
        disk_output[Band(1)] ==
        permutedims(permuted_output, (X, Y))
    @test abs(step(dims(raster_output, Y))) ≈
        abs(step(dims(raster_output, X))) ≈ 
        abs(step(dims(disk_output, X))) ≈ 
        abs(step(dims(permuted_output, X))) ≈ output_res
    @test name(cea) == name(raster_output)

    rm("resample.tif")

    @testset "missingval propagates" begin
        @test missingval(resample(cea; res=output_res, crs=output_crs, method)) == 0x00
    end

    @testset "snapped size and dim index match" begin
        snaptarget = raster_output
        snapped = resample(cea; to=snaptarget)
        disk_snapped = resample(cea; to=snaptarget, filename="raster.tif")
        @test size(snapped) == size(disk_snapped) == size(snaptarget)
        @test isapprox(index(snaptarget, Y), index(snapped, Y))
        @test isapprox(index(snaptarget, X), index(snapped, X))
        @test isapprox(index(snaptarget, Y), index(disk_snapped, Y))
        @test isapprox(index(snaptarget, X), index(disk_snapped, X))
    end

    @testset "`method` only does nothing" begin
        resampled = resample(cea; method)
        @test crs(cea) == crs(resampled)
        @test cea == resampled
        @test extent(cea) == extent(resampled)
    end

    @testset "only `res` kw changes the array size predictably" begin
        res = step(span(cea, X)) / 2
        resampled = resample(cea; res)
        @test crs(cea) == crs(resampled)
        @test size(dims(resampled, (X, Y))) == size(dims(cea, (X, Y))) .* 2
        # GDAL fp error see above
        # @test_broken extent(cea) == extent(resampled)
        resampled = resample(cea; res=(res, 2res))
        @test size(dims(resampled, (X, Y))) == (size(cea, X) .* 2, size(cea, Y))
        resampled = resample(cea; res=(X(2res), Y(res)))
        @test size(dims(resampled, (X, Y))) == (size(cea, X), size(cea, Y) * 2)
    end

    @testset "only size or res allowed not both" begin
        @test_throws ArgumentError resample(cea; res=output_res, size=(1000, 1000))
    end

    @testset "only `size` kw sets the size" begin
        res = step(span(cea, X)) / 2
        resampled = resample(cea; size=(100, 200))
        @test crs(cea) == crs(resampled)
        @test size(dims(resampled, (X, Y))) == size(resampled[:, :, 1]) == (100, 200)
        resampled = resample(cea; size=(X(99), Y(111)))
        @test crs(cea) == crs(resampled)
        @test size(dims(resampled, (X, Y))) == size(resampled[:, :, 1]) == (99, 111)
        resampled = resample(cea; size=100)
        @test size(dims(resampled, (X, Y))) == size(resampled[:, :, 1]) == (100, 100)
    end

    @testset "Extent `to` can resize arbitrarily" begin
        to = Extent(X=(-10000.0, 2000.0), Y=(4.2e6, 4.3e6))
        resampled = resample(cea; to)
        @test all(map((bs...,) -> all(map(≈, bs...)), extent(resampled), to))
        @test crs(cea) == crs(resampled)
        # Most of the area is extended, and filled with missingval
        @test sum(x -> x === missingval(resampled), resampled) > length(resampled) / 2
    end

    @testset "Geometries work as `to`" begin
        geom = GeoInterface.Polygon([[(0.0, 4e6), (-1e5, 4.4e6), (-1e5, 4.2e6)]])
        resampled = resample(cea; to=geom)
        @test map(extent(resampled[Band(1)]), GeoInterface.extent(geom)) do bs1, bs2 
            map(≈, bs2, bs2) |> all 
        end |> all
        @test crs(cea) == crs(resampled)
        # Most of the area is extended, and filled with missingval
        @test sum(x -> x === missingval(resampled), resampled) > length(resampled) / 2
    end

    @testset "only `crs` kw changes the array size" begin
        resampled = resample(cea; crs=EPSG(3857), method)
        @test size(dims(resampled, (X, Y))) !== size(dims(cea, (X, Y)))
        @test crs(resampled) == EPSG(3857)
    end

    @testset "no existing crs warns" begin
        nocrs = setcrs(cea, nothing)
        @test crs(nocrs) == nothing
        @test_warn "does not have crs" resample(nocrs; crs=output_crs, method)
    end

    @testset "resample eltype propagates" begin
        r = Raster(rand(UInt8, X(1:10), Y(1:10)))
        r1 = resample(r; to=r)
        @test eltype(r1) == UInt8
    end

    @testset "dimensions matcha after resampling with only `to`" begin
        # some weird dimensions
        to = (
            Lat(Projected(1:10, span = Regular(1 + eps()), crs = nothing, order = ForwardOrdered(), sampling = Intervals(Center()))),
            X(Sampled([1,2,3], span = Regular(1), order = ForwardOrdered(), sampling = Intervals(Start())))
        )

        resampled = resample(cea; to)
        @test dims(resampled) == to

        # test with 3d
        resampled_3D = resample(cat(cea, cea; dims = Z(1:2)); to)
        @test length(dims(resampled_3D)) == 3
        @test dims(resampled_3D, (1,2)) == to
        @test dims(resampled_3D, Z) == Z(1:2)
    end

    @testset "resample to the same size or resolution leaves raster unchanged" begin
        res = 2
        ys = range(; start=(90 - res / 2), step=-res, stop=(-90 + res / 2))
        xs = range(; start=(-180 + res / 2), step=res, stop=(180 - res / 2))

        point_dims = Y(ys), X(xs)
        interval_dims = Y(ys; sampling=Intervals(Center())), X(xs; sampling=Intervals(Center()))
        rev_y_point_dims = Y(reverse(ys)), X(xs)
        rev_y_interval_dims = reverse(interval_dims[1]), interval_dims[2]
        rev_x_point_dims = Y(ys), X(reverse(xs))
        rev_x_interval_dims = interval_dims[1], reverse(interval_dims[2])
        test_dims = (point_dims, interval_dims, rev_x_point_dims, rev_x_interval_dims, rev_y_point_dims, rev_y_interval_dims)
        ds_fwd = point_dims; f = identity
        ds_fwd = point_dims; f = reverse

        ds = ds_fwd
        raster = Raster(rand(ds), crs=EPSG(4326), missingval=NaN)
        resampled_res = resample(raster; res)
        resampled_size = resample(raster; size=size(raster), method=:near)
        @test resampled_res == resampled_size == raster
        @test dims(resampled_res) == dims(resampled_size) == dims(raster)

        ds = reverse(ds_fwd)
        raster = Raster(rand(ds), crs=EPSG(4326), missingval=NaN)
        resampled_res = resample(raster; res)
        resampled_size = resample(raster; size=size(raster), method=:near)
        @test resampled_res == resampled_size == raster
        @test dims(resampled_res) == dims(resampled_size) == dims(raster)
    end

end
