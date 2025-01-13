using Rasters, Test, Dates, DiskArrays
using Rasters.Lookups, Rasters.Dimensions
using Rasters: isdisk, ismem, filename

data1 = cumsum(cumsum(ones(10, 11); dims=1); dims=2)
data2 = 2cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = X(10:10:100), Y(-50:10:50) 
dims2 = (dims1..., Ti([DateTime(2019)]))
refdimz = ()
mval = -9999.0
meta = NoMetadata()
nme = :test

ga2 = Raster(data2, dims2)
ga1 = Raster(data1; dims=dims1, refdims=refdimz, name=nme, metadata=meta, missingval=mval)

@test ga1 == data1
@test ga2 == data2

@testset "from file" begin
    @test_throws ArgumentError Raster("notafile")
    @test_throws ArgumentError Raster("notafile", dims1)
end

@testset "Single dim vector constructor" begin
    A = [1, 2, 3]
    @test Raster(A, X()) === Raster(A, (X(),))
end

@testset "array properties" begin
    @test name(ga1) == :test
    @test missingval(ga1) == -9999.0
    @test isdisk(ga1) == false
    @test ismem(ga1) == true
    @test filename(ga1) === nothing
    @test DiskArrays.haschunks(ga1) == DiskArrays.Unchunked()
    @test DiskArrays.eachchunk(ga1) == DiskArrays.GridChunks((10, 11), (10, 11))
end

@testset "crs" begin
    crs(ga1) == nothing
    mappedcrs(ga1) == nothing
    gapr = setcrs(ga1, EPSG(4326))
    gampr = setmappedcrs(gapr, EPSG(4326))
    @test crs(gapr) == crs(gapr[X(1)]) == crs(dims(gapr)) == EPSG(4326)
    @test mappedcrs(gampr) == mappedcrs(gampr[X(1)]) == EPSG(4326)
end

@testset "array dims have been formatted" begin
    @test index(ga2) == (10.0:10:100, -50.0:10:50.0, [DateTime(2019)])
    @test dims(ga1)[1:2] == dims(ga2)[1:2]
end

@testset "constructors" begin
    x = X(Projected(1.0:1.0:2.0; crs=ProjString("+proj=")))
    y = Y(Projected(1.0:1.0:2.0; crs=ProjString("+proj=")))

    A1 = rand(Float32, x, y)
    @test A1 isa Raster{Float32,2}
    @test size(A1) == (2, 2)

    A2 = zeros(Int, x, y)
    @test A2 isa Raster{Int,2}
    @test A2 == [0 0; 0 0] 

    A3 = ones(x, y)
    @test A3 isa Raster{Float64,2}
    @test A3 == [1.0 1.0; 1.0 1.0] 

    A4 = fill(99.0, x, y)
    @test A4 isa Raster{Float64,2}
    @test A4 == [99.0 99.0; 99.0 99.0] 

    A5 = trues(x, y)
    @test A5 isa Raster{Bool,2}
    @test A5 == [true true; true true] 

    A6 = falses(x, y)
    @test A6 isa Raster{Bool,2}
    @test A6 == [false false; false false] 

    A7 = Raster(undef, extent(A1); size=(2, 2))
    @test A7 isa Raster{Float64,2}
    @test size(A7) == (2, 2)
    @test dims(A7) == dims(A1)

    A8 = Raster{Int}(undef, extent(A1); res=1.0)
    @test A8 isa Raster{Int,2}
    @test size(A8) == (2, 2)
    @test dims(A8) == dims(A1)

    A9 = Raster{Bool}(undef, extent(A1); res=1.0, sampling=Intervals(Center()), closed=false)
    @test A9 isa Raster{Bool,2}
    @test size(A9) == (2, 2)
    # Has to the same lookups to include the bounds in the extent
    @test collect.(dims(A9)) == collect.(dims(A1))
    @test sampling(A9, X) == Intervals(Center())

    A10 = Raster{Bool}(undef, extent(A1); size=(2, 2), sampling=Intervals(Center()))
    @test A10 isa Raster{Bool,2}
    @test size(A10) == (2, 2)
    # Has to the same lookups to include the bounds in the extent
    @test collect.(dims(A9)) == collect.(dims(A1))
    @test sampling(A10, X) == Intervals(Center())
end


@testset "collect and Array" begin
    @test collect(ga1) isa Array
    @test collect(ga1) == data1
    @test Array(ga1) isa Array{Float64,2}
    @test Array(ga1) == data1
end

@testset "show" begin
    sh = sprint(show, MIME("text/plain"), ga1)
    # Test but don't lock this down too much
    @test occursin("Raster", sh)
    @test occursin("Y", sh)
    @test occursin("X", sh)
end

@testset "skipmissing uses missingval" begin
    # Test missingval=NaN
    raster = Raster([NaN 1.0; 2.0 NaN], (X, Y); missingval=NaN)
    @test collect(skipmissing(raster)) == [2.0, 1.0]
    @test collect(keys(skipmissing(raster))) == [CartesianIndex(2, 1), CartesianIndex(1, 2)]
    @test collect(eachindex(skipmissing(raster))) == [2, 3]
    @test_throws MissingException skipmissing(raster)[1]
    @test skipmissing(raster)[2] == 2.0

    # It skips actual missing values as well
    mraster = Raster([NaN 1.0; missing NaN], (X, Y); missingval=NaN)
    @test collect(skipmissing(mraster)) == [1.0]
    @test collect(keys(skipmissing(mraster))) == [CartesianIndex(1, 2)]
    @test collect(eachindex(skipmissing(mraster))) == [3]
    @test skipmissing(mraster)[3] == 1.0
    @test_throws MissingException skipmissing(mraster)[2]

    # Test missingval=nothing
    fraster = Raster([NaN 1.0; 2.0 NaN], (X, Y); missingval=nothing)
    mraster = Raster([NaN 1.0; missing NaN], (X, Y); missingval=nothing)
    iraster = Raster([1 1; missing 2], (X, Y); missingval=nothing)
    nraster = Raster([1 1; missing nothing], (X, Y); missingval=nothing)
    @test length(collect(skipmissing(fraster))) == 4 # Keeps NaN
    @test length(collect(skipmissing(mraster))) == 3 # Drops missing
    @test length(collect(skipmissing(iraster))) == 3 # Drops missing (integers)
    @test length(collect(skipmissing(nraster))) == 3 # Keeps nothing (integer)

    # Confirm that missingvals are removed by value, even when types don't match
    r = Raster(ones(Int16, 8, 8), (X,Y); missingval=Int16(-9999))
    @test missingval(r) == Int16(-9999)
    r[1:4, 1:4] .= -9999
    rf = Float64.(r)
    @test missingval(rf) === -9999.0
    @test !(missingval(rf) in skipmissing(rf))
    @test length(collect(skipmissing(r))) == 48
end

@testset "table" begin
    table = DimTable(ga1)
    tra = Raster(table, dims(ga1))
    @test tra == ga1
    @test name(tra) == :test
    @test_throws ArgumentError Raster(table, dims(ga1); name=:x)
end

@testset "from vector" begin
    vra = Raster(vec(ga1), dims(ga1); name=:from_vector)
    @test vra == ga1
    @test name(vra) == :from_vector
end

@testset "keywords" begin
    md = Dict("test" => 1)
    rast = Raster(rand(X(10:0.1:20), Y(9:0.1:19)); 
        name=:test, missingval=9999.9, metadata=md,
        refdims=(Ti(DateTime(2000,1,1)),),
        crs=EPSG(4326), mappedcrs=EPSG(3857), 
    )
    @test name(rast) == :test
    @test missingval(rast) == 9999.9
    @test Rasters.metadata(rast) == md
    @test crs(rast) == EPSG(4326)
    @test mappedcrs(rast) == EPSG(3857)
    @test refdims(rast) == (Ti(DateTime(2000,1,1)),)
end
