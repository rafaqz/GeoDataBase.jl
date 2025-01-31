using Rasters, Test, Dates, Statistics
using Rasters.Lookups, Rasters.Dimensions
using Rasters: upsample, downsample, aggregate

@testset "upsample" begin
    @test upsample(1, 2) == 1
    @test upsample(2, 2) == 3
    @test upsample(3, 2) == 5
    @test upsample(1, 3) == 1
    @test upsample(2, 3) == 4
    @test upsample(3, 3) == 7
end

@testset "downsample" begin
    @test downsample(1, 2) == 1
    @test downsample(3, 2) == 2
    @test downsample(5, 2) == 3
    @test downsample(1, 3) == 1
    @test downsample(4, 3) == 2
    @test downsample(7, 3) == 3
end

data1 = [ 1  2  3  4  5  6 -1
          7  8  9 10 11 12 -1
         13 14 15 16 17 18 -1]
data2 = 2 * data1
data3 = 3 * data1
data4 = 4 * data1
dimz = X(Projected([30., 40., 50.]; order=ForwardOrdered(), span=Regular(10.0), sampling=Points(), crs=EPSG(4326))), 
       Y(Projected(LinRange(-10., 20., 7); order=ForwardOrdered(), span=Regular(5.0), sampling=Points(), crs=EPSG(4326)))
array1 = Raster(data1, dimz)
array2 = Raster(data2, dimz)
array1a = Raster(data3, dimz)
array2a = Raster(data4, dimz)
stack1 = RasterStack(array1, array2; name=(:array1, :array2))
stack2 = RasterStack(array1a, array2a; name=(:array1, :array2))
dates = DateTime(2017):Year(1):DateTime(2018)
series = RasterSeries([stack1, stack2], (Ti(dates),))

@testset "Aggregate a dimension" begin
    lat = Y(Sampled(LinRange(3, 13, 6), ForwardOrdered(), Regular(2.0), Intervals(Start()), NoMetadata()))
    aglat = aggregate(Start(), lat, 3)
    @test span(lookup(aglat)) == Regular(6.0)
    @test disaggregate(Start(), aglat, 3) == lat

    aglon = aggregate(Start(), dimz[1], 3)
    @test step(lookup(aglon)) === 30.0
    @test val(aglon) == [30.0]
    disaglon = disaggregate(Start(), aglon, 3)
    @test index(disaglon) == index(dimz[1])
    @test span(disaglon) == span(dimz[1])
    @test sampling(disaglon) == sampling(dimz[1])

    aglat = aggregate(Start(), dimz[2], 3)
    @test step(lookup(aglat)) === 15.0
    @test index(aglat) == LinRange(-10.0, 5.0, 2)
    disaglat = disaggregate(Start(), aglat, 3)
    # The last item is lost due to rounding in `aggregate`
    @test index(disaglat) != index(dimz[2])
    @test index(disaglat) === LinRange(-10.0, 15.0, 6)
    @test span(disaglat) == span(dimz[2])
    @test sampling(disaglat) == sampling(dimz[2])
end

@testset "aggregate a single dim" begin 
    ag = aggregate(Start(), series, (X(3),))
    @test size(first(ag)) == (1, 7)
    ag = aggregate(Start(), series, (Y(5),))
    @test size(first(ag)) == (3, 1)
    ag = aggregate(Start(), series, (Y(2),))
    @test size(first(ag)) == (3, 3)
end

@testset "dimension lookups are maintained" begin 
    x, y, z = X(1:100), Y(1:100), Z(["bar", "foo"])
    A = rand(x, y, z) |> Raster
    A_agg = Rasters.aggregate(sum, A, (X(5),))
    A_disagg = Rasters.disaggregate(A, (X(3),))
    @test dims(A_agg, (Y, Z)) == (y, z)
    @test dims(A_disagg, (Y, Z)) == (y, z)
end

@testset "aggregate and disaggregate at a locus" begin
    @testset "single scale single locus" begin
        scale = 3
        array1_aggstart = aggregate(Start(), array1, scale) 
        @test array1_aggstart == [1 4]
        @test length.(dims(array1_aggstart)) == size(array1_aggstart)
        array1_disagstart = disaggregate(Start(), array1_aggstart, scale) 
        disaggstart = 
            [1 1 1 4 4 4 
             1 1 1 4 4 4
             1 1 1 4 4 4]
        @test array1_disagstart == disaggstart
        array1_aggcenter = aggregate(Center(), array1, scale)
        @test array1_aggcenter == [8 11]
        array1_disagcenter = disaggregate(Center(), array1_aggcenter, scale)
        disaggcenter = 
            [8 8 8 11 11 11 
             8 8 8 11 11 11
             8 8 8 11 11 11]
        @test array1_disagcenter == disaggcenter
        array1_aggend = aggregate(End(), array1, scale) 
        @test array1_aggend == [15 18]
        array1_disaggend = disaggregate(End(), array1_aggend, scale)
        @test array1_disaggend == 
            [15 15 15 18 18 18 
             15 15 15 18 18 18
             15 15 15 18 18 18]
        stack1_aggstart = aggregate(Start(), stack1, scale)
        @test stack1_aggstart[:array2] == [2 8]
        stack1_disaggstart = disaggregate(Start(), stack1_aggstart, scale)
        @test stack1_disaggstart[:array1] == disaggstart 
        @test aggregate(Center(), stack1, scale)[:array2] == [16 22]
        series_aggcenter = aggregate(Center(), series, scale)
        @test series_aggcenter[2][:array2] == [32 44]
        series_disaggcenter = disaggregate(Center(), series_aggcenter, scale)
        @test series_disaggcenter[2][:array2] == 
            [32 32 32 44 44 44
             32 32 32 44 44 44
             32 32 32 44 44 44]
        @test typeof(aggregate(Start(), series, scale)) <: RasterSeries
    end

    @testset "mixed scales" begin
        scale = (3, 2)
        @test aggregate(Start(), array1, scale) == [1 3 5]
        scale = (Y(2), X(3))
        @test aggregate(Center(), array1, scale) == [8 10 12]
        @test aggregate(End(), array1, scale) == [14 16 18]
        A = aggregate(Start(), array1, scale)
        @test length.(dims(A)) == size(A)
    end

    @testset "mixed locus" begin
        scale = 3
        @test aggregate((End(), Start()), array1, 3) == [13 16]
        @test aggregate((End(), Start()), array1, (3, 2)) == [13 15 17]
        A = aggregate((End(), Start()), array1, scale)
        @test length.(dims(A)) == size(A)
    end

    @testset "dim scale" begin
        agg = aggregate(Start(), array1, (Y(3), X(1))) 
        @test agg == aggregate(Start(), array1, (1, 3))
        disagg = disaggregate(Start(), agg, (Y(3), X(1))) 
        @test disagg ==
            [1  1  1  4  4  4
             7  7  7 10 10 10
            13 13 13 16 16 16]
        agg = aggregate(Start(), array1, (X(1), Y(Near(-4)))) 
        @test agg == aggregate(Start(), array1, (1, 2))
        @testset "scale 1 dims are unchanged" begin
            @test dims(agg, X) == dims(array1, X)
        end
    end

end

@testset "Aggregate with Colon" begin
    @test aggregate(sum, array1, :) == [sum(array1);;]
    @test aggregate(sum, array1, (1, :)) == sum(array1; dims=2)
    @test aggregate(sum, array1, (X=1, Y=:)) == sum(array1; dims=2)
end

@testset "Aggregate with a function" begin
    @test aggregate(sum, array1, 3) == [72 99]
    @test aggregate(median, array1, 3) == [8 11]
    @test aggregate(sum, array1, (3, 2)) == [45 57 69]
    A = aggregate(sum, array1, (3, 2))
    @test length.(dims(A)) == size(A)
end

@testset "test all scale modes" begin
    @test aggregate(sum, array1, (3, 2)) == 
        aggregate(sum, array1, (Y(2), X(3))) == 
        aggregate(sum, array1, (Y=2, X=3)) == 
        aggregate(sum, array1, (Y=>2, X=>3)) == 
        aggregate(sum, array1, (:Y=>2, :X=>3)) == 
        [45 57 69]
end

@testset "Aggregate with a function with missing values" begin
    data_m = [ 1  2  3  4  5  6 -1
               7  8  9 10 11 12 -1
              13 14 15 16 missing 18 -1]
    array_m = Raster(data_m, dimz)
    @test all(aggregate(sum, array_m, 3) .=== [72 missing])
    @test all(aggregate(sum, array_m, 3; skipmissing=true) .=== [72 82])
    data_m0 = [ 1  2  3  4  5  6 -1
                7  8  0 10 11 12 -1
               13 14 15 16 17 18 -1]
    array_m0 = Raster(data_m0, dimz; missingval=0)
    @test aggregate(sum, array_m0, 3) == [0 99]
    @test aggregate(sum, array_m0, 3; skipmissing=true) == [63 99]
    @test all(aggregate(mean, array_m0, 3) .=== [0.0 11.0])
end

@testset "Aggregate to an array with a different missing value" begin
    data_m = [ 1  2  3  4  5  6 -1
              7  8  9 10 11 12 -1
             13 14 15 16 missing 18 -1]
    src = Raster(data_m, dimz)
    dst = Raster(zeros(1, 2), map(d -> aggregate(mean, d, 3), dimz); missingval=-9999.0)
    @test all(aggregate!(sum, dst, src, 3) .=== [72.0 -9999.0])
end

@testset "Aggregate different index lookups" begin
    dimz = Y([1, 3, 2]), Dim{:category}([:a, :b, :c]), X([10, 20, 30, 40])
    a1 = [1 2 3; 4 5 6; 7 8 9]
    A = cat(a1, a1 .+ 10, a1 .+ 20, a1 .+ 30, dims=3)
    da = Raster(A, dimz)
    @test vec(aggregate(sum, da, (3, 2, 2))) == [114, 354]
    @test size(aggregate(sum, da, 3)) == (3, 3, 1) 
end

@testset "Aggregate ignores categorical by default" begin
    rast = Raster(rand(X(1:10), Y(1:10), Z([:a, :b, :c]))) 
    skipped_ag = @test_logs (:info,) aggregate(sum, rast, 2)
    @test size(skipped_ag) == (5, 5, 3)
    # Unless specified explicitly
    explicit_ag = @test_nowarn aggregate(sum, rast, (X=5, Y=5, Z=3))
    @test size(explicit_ag) == (2, 2, 1)
    # Aggregated Categorical become NoLookup
    @test isnolookup(explicit_ag, Z)
    @test_throws ArgumentError aggregate(sum, rast[X=1, Y=1], 2)
end

