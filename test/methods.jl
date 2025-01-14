using Rasters, Test, ArchGDAL, ArchGDAL.GDAL, Dates, Statistics, DataFrames, Extents, Shapefile, GeometryBasics
import GeoInterface
using Rasters.Lookups, Rasters.Dimensions 
using Rasters: bounds, trim

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

A = [missing 7.0f0; 2.0f0 missing]
B = [1.0 0.4; 2.0 missing]
ga = Raster(A, (X(1.0:1:2.0), Y(1.0:1:2.0)); missingval=missing) 
st = RasterStack((a=A, b=B), (X, Y); missingval=(a=missing,b=missing))
st2 = RasterStack((a=A[1,:], b=B), (X, Y); missingval=(a=missing,b=missing))
se = RasterSeries([ga, ga], Rasters.Band(1:2))

pointvec = [(-20.0, 30.0),
            (-20.0, 10.0),
            (0.0, 10.0),
            (0.0, 30.0),
            (-20.0, 30.0)]

            
vals = [1, 2, 3, 4, 5]
polygon = ArchGDAL.createpolygon(pointvec)
multi_polygon = ArchGDAL.createmultipolygon([[pointvec]])
multi_polygon = ArchGDAL.createmultipolygon([[pointvec]])
multi_point = ArchGDAL.createmultipoint(pointvec)
linestring = ArchGDAL.createlinestring(pointvec)
multi_linestring = ArchGDAL.createmultilinestring([pointvec])
linearring = ArchGDAL.createlinearring(pointvec)
pointfc = map(GeoInterface.getpoint(polygon), vals) do geom, v
    (geometry=geom, val1=v, val2=2.0f0v)
end
polytemplate = Raster(ones(X(-20:5; sampling=Intervals(Center())), Y(0:30; sampling=Intervals(Center()))))

test_shape_dir = realpath(joinpath(dirname(pathof(Shapefile)), "..", "test", "shapelib_testcases"))
shp_paths = filter(x -> occursin("shp", x), readdir(test_shape_dir; join=true))
shppath = shp_paths[1]
shphandle = Shapefile.Handle(shppath)

ga99 = replace_missing(ga, -9999)
gaNaN = replace_missing(ga, NaN32)
gaMi = replace_missing(ga)

@testset "replace_missing" begin
    @test eltype(ga99) == Float32
    @test eltype(gaNaN) == Float32
    @test eltype(gaMi) == Union{Float32,Missing}
    @test eltype(replace_missing(ga, 0.0)) == Float64
    @test all(isequal.(ga99, [-9999.0f0 7.0f0; 2.0f0 -9999.0f0]))
    @test missingval(ga99) === -9999.0f0
    @test all(isequal.(gaNaN, [NaN32 7.0f0; 2.0f0 NaN32]))
    @test missingval(gaNaN) === NaN32
    @test all(isequal.(gaMi, ga))
    @test missingval(gaMi) === missing
    # The second layer NaN32 will be promoted to NaN to match the array type
    @test missingval(replace_missing(st, NaN32)) === (a=NaN32, b=NaN)
    @test all(map(values(replace_missing(st, NaN32)), (a=[NaN32 7.0f0; 2.0f0 NaN32], b=[1.0 0.4; 2.0 NaN])) do x, y
        all(x .=== y)
    end)
    testfile = tempname() * ".tif"
    dNaN = replace_missing(ga, NaN32; filename=testfile)
    read(dNaN)
    @test all(isequal.(dNaN, [NaN32 7.0f0; 2.0f0 NaN32]))
    stNaN = replace_missing(st, NaN32; filename="teststack.tif")
    @test all(maplayers(stNaN[Band(1)], (a=[NaN32 7.0f0; 2.0f0 NaN32], b=[1.0 0.4; 2.0 NaN])) do x, y
        all(x .=== y)
    end)
    rm("teststack_a.tif")
    rm("teststack_b.tif")
end

@testset "boolmask" begin
    @test boolmask(ga) == [false true; true false]
    @test boolmask(ga; invert=true) == [true false; false true]
    @test parent(boolmask(ga)) isa BitMatrix
    @test boolmask(ga99) == [false true; true false]
    @test boolmask(ga99; invert=true) == [true false; false true]
    @test boolmask(gaNaN) == [false true; true false]
    @test boolmask(gaNaN; invert=true) == [true false; false true]
    @test all(boolmask(st[(:b, :a)], alllayers=true) .=== [false true; true false])
    @test all(boolmask(st[(:b, :a)], alllayers=true, invert=true) .=== [true false; false true])
    @test all(boolmask(st[(:b, :a)], alllayers=false) .=== [true true; true false])    
    @test all(boolmask(st[(:b, :a)], alllayers=false, invert=true) .=== [false false; false true])    
    @test all(boolmask(st[(:b, :a)], alllayers=false, missingval=7.0) .=== [true true; true true])    
    @test all(boolmask(st[(:b, :a)], alllayers=false, missingval=7.0, invert=true) .=== [false false; false false])    
    @test all(boolmask(st[(:b, :a)], alllayers=true, missingval=7.0) .=== [true false; true true])    
    @test all(boolmask(st[(:b, :a)], alllayers=true, missingval=7.0, invert=true) .=== [false true; false false])    
    @test all(boolmask(st, alllayers=true, missingval=(a=missing, b=0.4)) .=== [false false; true false])    
    @test_throws ArgumentError boolmask(st, alllayers=true, missingval=(b=missing, a=0.4))  
    se2 = RasterSeries([st.b, st.a], Rasters.Band(1:2))
    @test all(boolmask(se2, alllayers=true) .=== [false true; true false])
    @test all(boolmask(se2, alllayers=false) .=== [true true; true false])    
    @test dims(boolmask(ga)) === dims(ga)
    x = boolmask(polygon; res=1.0, boundary=:touches) 
    @test x == trues(X(Projected(-20:1.0:0.0; sampling=Intervals(Start()), crs=nothing)),
                     Y(Projected(10.0:1.0:30.0; sampling=Intervals(Start()), crs=nothing)))
    @test all(x .!= boolmask(polygon; res=1.0, invert=true, boundary=:touches))
    @test parent(x) isa BitMatrix
    # With a :geometry axis
    x1 = boolmask([polygon, polygon]; collapse=false, res=1.0, boundary=:touches)
    x2 = boolmask([polygon, polygon]; collapse=false, res=1.0, boundary=:touches, sampling=Intervals(Center()))
    x3 = boolmask([polygon, polygon]; collapse=false, res=1.0, boundary=:touches, sampling=Intervals(End()))
    @test extent(x1) == extent(x2) == extent(x3) == Extent(X = (-20.0, 1.0), Y = (10.0, 31.0), geometry = (1, 2))
    x4 = boolmask([polygon, polygon]; collapse=false, size=(21, 21), boundary=:touches)
    x5 = boolmask([polygon, polygon]; collapse=false, size=(21, 21), boundary=:touches, sampling=Intervals(Center()))
    x6 = boolmask([polygon, polygon]; collapse=false, size=(21, 21), boundary=:touches, sampling=Intervals(End()))
    xs = (x1, x2, x3, x4, x5, x6)
    @test all(x1 .!= boolmask([polygon, polygon]; collapse=false, res=1.0, invert=true, boundary=:touches))
    @test sampling(x1, X) isa Intervals{Start}
    @test sampling(x2, X) isa Intervals{Center}
    @test sampling(x3, X) isa Intervals{End}
    @test sampling(x4, X) isa Intervals{Start}
    @test sampling(x5, X) isa Intervals{Center}
    @test sampling(x6, X) isa Intervals{End}
    for x in xs
        @test eltype(x1) == Bool
        @test size(x) == (21, 21, 2)
        @test sum(x) == 882
        @test parent(x) isa BitArray{3}
        @test eltype(x) == Bool
    end
    @testset "size adds nextfloat" begin
        s = boolmask([polygon, polygon]; collapse=false, size=(21, 21), boundary=:touches)
        bounds(s, X)[1] == -20.0
        bounds(s, X)[2] > 0.0
        bounds(s, Y)[1] == 10.0
        bounds(s, Y)[2] > 30.0
        c = boolmask([polygon, polygon]; collapse=false, size=(21, 21), boundary=:touches, sampling=Intervals(Center()))
        bounds(c, X)[1] == -20.0
        bounds(c, X)[2] > 0.0
        bounds(c, Y)[1] == 10.0
        bounds(c, Y)[2] > 30.0
        e = boolmask([polygon, polygon]; collapse=false, size=(21, 21), boundary=:touches, sampling=Intervals(End()))
        bounds(e, X)[1] < -20.0
        bounds(e, X)[2] == 0.0
        bounds(e, Y)[1] < 10.0
        bounds(e, Y)[2] == 30.0
    end

    x = boolmask([polygon, polygon]; collapse=true, res=1.0, boundary=:touches)
    @test all(x .!= boolmask([polygon, polygon]; collapse=true, res=1.0, invert=true, boundary=:touches))
    @test size(x) == (21, 21)
    @test sum(x) == 441
    @test parent(x) isa BitMatrix
    for poly in (polygon, multi_polygon) 
        @test boolmask(poly; to=polytemplate) == .!boolmask(poly; to=polytemplate, invert=true)
        @test boolmask(poly; to=polytemplate, shape=:line) == .!boolmask(poly; to=polytemplate, shape=:line, invert=true)
        @test boolmask(poly; to=polytemplate, shape=:point) == .!boolmask(poly; to=polytemplate, shape=:point, invert=true)
    end
    # TODO: use explicit intervals in Extents.jl to make this exact?
    @testset "slightly larger extent" begin
        rast = boolmask([polygon, polygon]; shape=:line, collapse=false, res=1.0)
        @test GeoInterface.extent(rast).X[1] == GeoInterface.extent(polygon).X[1]
        @test GeoInterface.extent(rast).X[2] > GeoInterface.extent(polygon).X[2]
        @test GeoInterface.extent(rast).Y[1] == GeoInterface.extent(polygon).Y[1]
        @test GeoInterface.extent(rast).Y[2] > GeoInterface.extent(polygon).Y[2]
    end
end

@testset "missingmask" begin
    @test all(missingmask(ga) .=== [missing true; true missing])
    @test all(missingmask(ga; invert=true) .=== [true missing; missing true])
    @test all(missingmask(ga99) .=== [missing true; true missing])
    @test all(missingmask(ga99; invert=true) .=== [true missing; missing true])
    @test all(missingmask(gaNaN) .=== [missing true; true missing])
    @test all(missingmask(gaNaN; invert=true) .=== [true missing; missing true])

    @test all(missingmask(st[(:b, :a)], alllayers=true) .=== [missing true; true missing])
    @test all(missingmask(st[(:b, :a)], alllayers=true, invert=true) .=== [true missing; missing true])
    @test all(missingmask(st[(:b, :a)], alllayers=false) .=== [true true; true missing])  
    @test all(missingmask(st[(:b, :a)], alllayers=false, missingval = 7.0) .=== [true true; true true])    
    @test all(missingmask(st[(:b, :a)], alllayers=true, missingval = 7.0) .=== [true missing; true true])      
    @test dims(missingmask(ga)) == dims(ga)
    @test all(missingmask(st[(:b, :a)], alllayers = true) .=== [missing true; true missing])
    @test all(missingmask(st[(:b, :a)], alllayers = false) .=== [true true; true missing])    
    mm_st2 = missingmask(st2)
    mm_st2_inverted = missingmask(st2; invert=true)
    @test dims(mm_st2) == dims(mm_st2_inverted) == dims(st2)
    @test all(mm_st2 .=== [missing missing; true missing])    
    @test all(mm_st2_inverted .=== [true true; missing true])    
    @test all(missingmask(st2, alllayers = false) .=== [missing; true])    
    @test all(missingmask(se) .=== missingmask(ga))
    @test missingmask(polygon; res=1.0, boundary=:touches) == 
        fill!(Raster{Union{Missing,Bool}}(undef, X(Projected(-20:1.0:0.0; crs=nothing)), Y(Projected(10.0:1.0:30.0; crs=nothing))), true)
    x = missingmask([polygon, polygon]; collapse=false, res=1.0, boundary=:touches)
    x_inverted = missingmask([polygon, polygon]; collapse=false, res=1.0, invert=true, boundary=:touches)
    @test all(ismissing.(x_inverted))
    @test eltype(x) == Union{Bool,Missing}
    @test size(x) == (21, 21, 2)
    @test sum(x) == 882
    @test parent(x) isa Array{Union{Missing,Bool},3}
    x = missingmask([polygon, polygon]; collapse=true, res=1.0, boundary=:touches)
    @test size(x) == (21, 21)
    @test sum(x) == 441
    @test parent(x) isa Array{Union{Missing,Bool},2}
    for poly in (polygon, multi_polygon) 
        @test all(missingmask(poly; to=polytemplate) .=== 
                  replace(missingmask(poly; to=polytemplate, invert=true), missing=>true, true=>missing))
        @test all(missingmask(poly; to=polytemplate, shape=:line) .=== 
                  replace(missingmask(poly; to=polytemplate, shape=:line, invert=true), missing=>true, true=>missing))
        @test all(missingmask(poly; to=polytemplate, shape=:point) .=== 
                  replace(missingmask(poly; to=polytemplate, shape=:point, invert=true), missing=>true, true=>missing))
    end
end

@testset "mask" begin
    A1 = [missing 1; 2 3]
    A2 = view([0 missing 1; 0 2 3], :, 2:3)
    ga1 = Raster(A1, (X, Y); missingval=missing)
    ga2 = Raster(A2, (X, Y); missingval=missing)
    @test all(mask(ga1; with=ga) .=== mask(ga2; with=ga) .=== [missing 1; 2 missing])
    @test all(mask(ga1; with=ga, invert=true) .=== mask(ga2; with=ga, invert=true) .=== [missing missing; missing 3])
    ga2 = replace_missing(ga1 .* 1.0; missingval=NaN)
    @test all(mask(ga2; with=ga) .=== [NaN 1.0; 2.0 NaN])
    @test all(mask(ga2; with=ga, invert=true) .=== [NaN NaN; NaN 3.0])
    ga3 = replace_missing(ga1; missingval=-9999)
    mask!(ga3; with=ga)
    @test all(ga3 .=== [-9999 1; 2 -9999])
    ga4 = replace_missing(ga1; missingval=-9999)
    mask!(ga4; with=ga, invert=true)
    @test all(ga4 .=== [-9999 -9999; -9999 3])
    maskfile = tempname() * ".tif"
    dmask = mask(ga3; with=ga, filename=maskfile)
    @test Rasters.isdisk(dmask)
    stmask = mask(replace_missing(st, NaN); with=ga, filename="mask.tif")
    @test Rasters.isdisk(stmask)
    rm("mask_a.tif")
    rm("mask_b.tif")
    poly = polygon
    @testset "to polygon" begin
        for poly in (polygon, multi_polygon) 
            st1 = RasterStack(polytemplate, polytemplate)
            ser1 = RasterSeries([polytemplate, polytemplate], Ti(1:2))
            @test all(
                mask(polytemplate; with=polygon) .===
                mask(st1; with=polygon)[:layer1] .===
                mask(ser1; with=polygon)[1]
            )
            @test all(
                mask(polytemplate; with=polygon, invert=true) .===
                mask(st1; with=polygon, invert=true)[:layer1] .===
                mask(ser1; with=polygon, invert=true)[1] .===
                replace(replace(mask(polytemplate; with=polygon), missing => 0.0, 1.0 => missing), 0.0 => 1.0)
            )
            @test all(
                mask(polytemplate; with=polygon, invert=true, shape=:line) .===
                mask(st1; with=polygon, invert=true, shape=:line)[:layer1] .===
                mask(ser1; with=polygon, invert=true, shape=:line)[1] .===
                replace(replace(mask(polytemplate; with=polygon, shape=:line), missing => 0.0, 1.0 => missing), 0.0 => 1.0)
            )
            # TODO: investigate this more for Points/Intervals
            # Exactly how do we define when boundary values are inside/outside a polygon
            @test sum(skipmissing(mask(polytemplate; with=polygon, boundary=:inside))) == 19 * 19
            @test sum(skipmissing(mask(polytemplate; with=polygon, boundary=:inside, invert=true))) == prod(size(polytemplate)) - 19 * 19
            @test sum(skipmissing(mask(polytemplate; with=polygon, boundary=:center))) == 20 * 20
            @test sum(skipmissing(mask(polytemplate; with=polygon, boundary=:center, invert=true))) == prod(size(polytemplate)) - 20 * 20
            @test sum(skipmissing(mask(polytemplate; with=polygon, boundary=:touches))) == 21 * 21
            @test sum(skipmissing(mask(polytemplate; with=polygon, boundary=:touches, invert=true))) == prod(size(polytemplate)) - 21 * 21
        end
    end

    @testset "geometry encompassing raster" begin
        geom = GeoInterface.Polygon([GeoInterface.LinearRing([(0.0, 0.0), (0.0, 10.0), (10.0, 10.0), (10.0, 0.0), (0.0, 0.0)])])
        raster = Raster(ones(X(1:0.1:2), Y(1:0.1:2)), missingval=false)
        @test sum(mask(raster; with=geom)) == sum(raster)
    end
end

@testset "mask" begin
    # Floating point rasters
    a = Raster([1.0 0.0; 1.0 1.0], dims=(X, Y), missingval=0.0)
    b = Raster([1.0 1.0; 1.0 0.0], dims=(X, Y), missingval=0.0)

    # Integer rasters
    c = Raster([1 0; 1 1], dims=(X, Y), missingval=0)
    d = Raster([1 1; 1 0], dims=(X, Y), missingval=0)

    # Test that missingval is replaced in source mask (Floats)
    @test isequal(mask(a; with=b, missingval=3.14), [1.0 3.14; 1.0 3.14]) # Test missingval = 3.14
    @test isequal(mask(a; with=b, missingval=3.14, invert=true), [3.14 3.14; 3.14 1.0]) # Test missingval = 3.14
    @test isequal(mask(a; with=b, missingval=missing), [1.0 missing; 1.0 missing]) # Test missingval = missing
    @test isequal(mask(a; with=b, missingval=missing, invert=true), [missing missing; missing 1.0]) # Test missingval = missing
    @test isequal(mask(a; with=b, missingval=NaN), [1.0 NaN; 1.0 NaN]) # Test missingval = NaN
    @test isequal(mask(a; with=b, missingval=NaN, invert=true), [NaN NaN; NaN 1.0]) # Test missingval = NaN
    @test isequal(mask(a; with=b, missingval=NaN32), [1.0 NaN; 1.0 NaN]) # Test convert NaN32 to NaN
    @test isequal(mask(a; with=b, missingval=Inf), [1.0 Inf; 1.0 Inf]) # Test missingval = Inf
    @test_throws MethodError mask(a, with=b, missingval=nothing)

    # Test that missingval is replaced in source mask (Ints)
    @test isequal(mask(c, with=d, missingval=missing), [1 missing; 1 missing]) # Test missingval = missing
    @test isequal(mask(c, with=d, missingval=missing, invert=true), [missing missing; missing 1]) # Test missingval = missing
    @test isequal(mask(c, with=d, missingval=-1.0), [1 -1; 1 -1])
    @test isequal(mask(c, with=d, missingval=-1.0), [1 -1; 1 -1])
    @test_throws MethodError mask(c, with=d, missingval=nothing)
    @test_throws InexactError mask(c, with=d, missingval=NaN)
    @test_throws InexactError mask(c, with=d, missingval=3.14)
    @test_throws InexactError mask(c, with=d, missingval=Inf)

    # Test Type Stability
    @test eltype(mask(a, with=b, missingval=0)) == Float64
    @test eltype(mask(a, with=b, missingval=-1)) == Float64
    @test eltype(mask(a, with=b, missingval=Inf32)) == Float64
    @test eltype(mask(Float32.(a), with=b, missingval=Inf)) == Float32
    @test eltype(mask(Float32.(a), with=b, missingval=NaN)) == Float32
    @test eltype(mask(Float32.(a), with=b, missingval=0.0)) == Float32
    @test eltype(mask(Float32.(a), with=b, missingval=0)) == Float32
    @test eltype(mask(Float32.(a), with=b, missingval=-1)) == Float32
    @test eltype(mask(c, with=d, missingval=-1.0)) == Int64
    @test eltype(mask(c, with=d, missingval=0.0f0)) == Int64
    @test eltype(mask(c, with=Float64.(d), missingval=-1.0)) == Int64
    @test eltype(mask(c, with=Float64.(d), missingval=0.0f0)) == Int64

    # Test mask!
    @test_throws MethodError mask!(a, with=b, missingval=missing)
    @test isequal(mask!(deepcopy(a), with=b, missingval=3.14), [1.0 3.14; 1.0 3.14]) # Test missingval = 3.14
    @test isequal(mask!(deepcopy(a), with=b, missingval=NaN), [1.0 NaN; 1.0 NaN]) # Test missingval = NaN
    @test isequal(mask!(deepcopy(a), with=b, missingval=NaN32), [1.0 NaN; 1.0 NaN]) # Test convert NaN32 to NaN
    @test isequal(mask!(deepcopy(a), with=b, missingval=Inf), [1.0 Inf; 1.0 Inf]) # Test missingval = Inf
    @test isequal(mask(deepcopy(c), with=d, missingval=-1.0), [1 -1; 1 -1])
    @test_throws MethodError mask!(deepcopy(a), with=b, missingval=missing)
    @test_throws ArgumentError mask!(deepcopy(c), with=d, missingval=nothing)
    @test_throws InexactError mask!(deepcopy(c), with=d, missingval=NaN)
    @test_throws InexactError mask!(deepcopy(c), with=d, missingval=3.14)
    @test_throws InexactError mask!(deepcopy(c), with=d, missingval=Inf)
end

@testset "zonal" begin
    a = Raster((1:26) * (1:31)', (X(-20:5), Y(0:30)))
    pointvec_empty = [(-100.0, 0.0), (-100.0, 0.0), (-100.0, 0.0), (-100.0, 0.0), (-100.0, 0.0)]
    polygon_empty = ArchGDAL.createpolygon(pointvec_empty)
    @test zonal(sum, a; of=polygon) ==
        zonal(sum, a; of=[polygon, polygon])[1] ==
        zonal(sum, a; of=[polygon, polygon_empty])[1] ==
        zonal(sum, a; of=(geometry=polygon, x=:a, y=:b)) ==
        zonal(sum, a; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
        zonal(sum, a; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
        sum(skipmissing(mask(a; with=polygon)))
    @test zonal(sum, a; of=a) == 
        zonal(sum, a; of=dims(a)) ==
        zonal(sum, a; of=Extents.extent(a)) == 
        sum(a)

    b = a .* 2
    c = cat(a .* 3, a .* 3; dims=:newdim)
    st = RasterStack((; a, b, c))
    @test zonal(sum, st; of=polygon) == zonal(sum, st; of=[polygon])[1] ==
        zonal(sum, st; of=(geometry=polygon, x=:a, y=:b)) ==
        zonal(sum, st; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
        zonal(sum, st; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
        maplayers(sum ∘ skipmissing, mask(st; with=polygon))  
    @test zonal(sum, st; of=st) == 
        zonal(sum, st; of=dims(st)) == 
        zonal(sum, st; of=Extents.extent(st)) == 
        sum(st)

    @testset "skipmissing" begin
        a = Array{Union{Missing,Int}}(undef, 26, 31)
        a .= (1:26) * (1:31)'
        a[1:10, 3:10] .= missing
        rast = Raster(a, (X(-20:5), Y(0:30)))
        @test zonal(sum, rast; of=polygon, skipmissing=false) === missing
        @test zonal(sum, rast; of=polygon, skipmissing=true) isa Int
        @test !zonal(x -> x isa Raster, rast; of=polygon, skipmissing=true)
        @test zonal(x -> x isa Raster, rast; of=polygon, skipmissing=false)
    end
end

@testset "zonal return missing" begin
    a = Raster((1:26) * (1:31)', (X(-20:5), Y(0:30)))
    out_bounds_pointvec = [(-40.0, -40.0), (-40.0, -35.0), (-35.0, -35.0), (-35.0, -40.0)]
    out_bounds_polygon = ArchGDAL.createpolygon(out_bounds_pointvec)
    @test ismissing(zonal(sum, a; of=[polygon, out_bounds_polygon, polygon])[2]) &&
        ismissing(zonal(sum, a; of=[out_bounds_polygon, polygon])[1]) &&
        ismissing(zonal(sum, a; of=(geometry=out_bounds_polygon, x=:a, y=:b))) &&
        ismissing(zonal(sum, a; of=[(geometry=out_bounds_polygon, x=:a, y=:b)])[1])
    @test zonal(sum, a; of=[out_bounds_polygon, out_bounds_polygon, polygon])[3] == 
        sum(skipmissing(mask(a; with=polygon)))
end

@testset "classify" begin
    A1 = [missing 1; 2 3]
    ga1 = Raster(A1, (X, Y); missingval=missing)
    @test all(classify(ga1, 1=>99, 2=>88, 3=>77) .=== [missing 99; 88 77])
    @test all(classify(ga1, 1=>99, 2=>88, 3=>77; others=0) .=== [missing 99; 88 77])
    @test all(classify(ga1, 1=>99, 2=>88; others=0) .=== [missing 99; 88 0])
    A2 = [1.0 2.5; 3.0 4.0]
    ga2 = Raster(A2 , (X, Y))
    @test classify(ga2, (2, 3)=>:x, >(3)=>:y) == [1.0 :x; 3.0 :y]
    @test classify(ga2, (>=(1), <(2))=>:x, >=(3)=>:y) == [:x 2.5; :y :y]
    classify!(ga2, (1, 2.5)=>0.0, >=(3)=>-1.0; lower=(>), upper=(<=))
    @test ga2 == [1.0 0.0; -1.0 -1.0]
    classify!(ga2, [1 2.5 0.0; 2.5 4.0 -1.0]; lower=(>), upper=(<=))
    @test ga2 == [1.0 0.0; -1.0 -1.0]
    @test all(classify(ga1, [1 99; 2 88; 3 77]) .=== [missing 99; 88 77])
    @test_throws ArgumentError classify(ga1, [1, 2, 3])
end

@testset "points" begin    
    dimz = (X(9.0:1.0:10.0), Y(0.1:0.1:0.2))
    rast = Raster([1 2; 3 4], dimz; name=:test)
    rast2 = Raster([5 6; 7 8], dimz; name=:test2, missingval=5)
    rast_m = Raster([1 2; 3 missing], dimz; name=:test)
    table = (geometry=[missing, (9.0, 0.1), (9.0, 0.2), (10.0, 0.3)], foo=zeros(4))
    st = RasterStack(rast, rast2)
    ga = Raster(A, (X(9.0:1.0:10.0), Y(0.1:0.1:0.2)))
    @test all(collect(points(ga; order=(Y, X))) .=== [missing (0.2, 9.0); (0.1, 10.0) missing])
    @test all(collect(points(ga; order=(X, Y))) .=== [missing (9.0, 0.2); (10.0, 0.1) missing])
    @test all(points(ga; order=(X, Y), ignore_missing=true) .===
              [(9.0, 0.1) (9.0, 0.2); (10.0, 0.1) (10.0, 0.2)])
end

createpoint(args...) = ArchGDAL.createpoint(args...)

@testset "trim, crop, extend" begin
    A = [missing missing missing
         missing 2.0     0.5
         missing 1.0     missing]

    r_fwd = Raster(A, (X(1.0:1.0:3.0), Y(1.0:1.0:3.0)); missingval=missing)
    r_crs = Raster(A, (X(1.0:1.0:3.0), Y(1.0:1.0:3.0)); missingval=missing, crs = EPSG(4326))
    r_mcrs = Raster(A, (X(1.0:1.0:3.0), Y(1.0:1.0:3.0)); missingval=missing, crs = EPSG(4326), mappedcrs = EPSG(4326))
    r_revX = reverse(r_fwd; dims=X)
    r_revY = reverse(r_fwd; dims=Y)
    st1 = RasterStack((a = r_fwd, b = r_fwd))
    st2 = RasterStack((a = r_fwd, b = r_fwd), crs = EPSG(4326))
    st3 = RasterStack((a = r_fwd, b = r_fwd), crs = EPSG(4326), mappedcrs = EPSG(4326))

    ga = r_fwd
    for ga in (r_fwd, r_crs, r_mcrs, r_revX, r_revY, st1, st2, st3)
        # Test with missing on all sides
        ga_r = rot180(ga)
        trimmed = trim(ga)
        trimmed_r = trim(ga_r)
        trimmed_d = trim(ga_r)
        if ga === r_fwd
            @test all(trimmed .=== [2.0 0.5; 1.0 missing])
            @test all(trimmed_r .=== [missing 1.0; 0.5 2.0])
        end
        cropped = crop(ga; to=trimmed)
        cropped1 = crop(crop(ga; to=dims(trimmed, X)); to=dims(trimmed, Y))
        _, cropped2 = crop(trimmed, ga)
        cropped_r = crop(ga_r; to=trimmed_r)
        @test map(identicalelements, maybe_layers.((cropped, cropped1, trimmed))...) |> all
        @test map(identicalelements, maybe_layers.((cropped_r, trimmed_r))...) |> all
        extended = extend(cropped, ga)[1]
        extended_r = extend(cropped_r; to=ga_r)
        extended1 = extend(extend(cropped; to=dims(ga, X)); to=dims(ga, Y))
        filename = tempname() * ".tif"
        extended_d = extend(cropped; to=ga, filename)
        @test map(identicalelements, maybe_layers.((extended, extended1, replace_missing(extended_d), ga))...) |> all
        @test map(identicalelements, maybe_layers.((extended_r, ga_r))...) |> all
        @test all(map(==, lookup(extended_d), lookup(extended)))
    end

    @testset "unformatted dimension works in crop" begin
        xdim, ydim = X(1.0:0.2:2.0), Y(1.0:1:2.0)
        A = Raster(rand((X(1.0:0.2:4.0), ydim)))
        @test lookup(crop(A; to=xdim), X) == 1.0:0.2:2.0

        A = Raster(rand((X(1.0:0.2:4.0), ydim)))
        @test lookup(crop(A; to=xdim), X) == 1.0:0.2:2.0

        A = Raster(ones((X(1.0:0.2:1.4), ydim)); missingval=0.0)
        extend(A; to=xdim)
        @test lookup(extend(A; to=xdim), X) == 1.0:0.2:2.0
        @test extend(A; to=xdim) == [1.0 1.0; 1.0 1.0; 1.0 1.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]  
    end

    @testset "to polygons" begin
        A1 = Raster(zeros(X(-20:-5; sampling=Points()), Y(0:30; sampling=Points())))
        A2 = Raster(ones(X(-20:-5; sampling=Points()), Y(0:30; sampling=Points())))
        A1crop1 = crop(A1; to=polygon)
        A1crop2, A2crop = crop(A1, A2; to=polygon)
        size(A1crop1)
        @test size(A1crop1) == size(A1crop2) == size(A2crop) == (16, 21)
        @test bounds(A1crop1) == bounds(A1crop2) == bounds(A2crop) == ((-20, -5), (10, 30))
        A1extend1 = extend(A1; to=polygon)
        A1extend2, A2extend = extend(A1, A2; to=polygon)
        @test all(A1extend1 .=== A1extend2)
        @test size(A1extend1) == size(A1extend2) == size(A2extend) == (21, 31)
        @test bounds(A1extend1) == bounds(A1extend2) == bounds(A2extend) == ((-20.0, 0.0), (0, 30))
    end
    @testset "to featurecollection and table" begin
        A1 = Raster(zeros(X(-20:-5; sampling=Points()), Y(0:30; sampling=Points())))
        featurecollection = map(GeoInterface.getpoint(polygon), vals) do geometry, v
            (; geometry, val1=v, val2=2.0f0v)
        end
        fccrop = crop(A1; to=featurecollection)
        table = DataFrame(featurecollection)
        tablecrop = crop(A1; to=table)
        @test size(fccrop) == size(tablecrop) == (16, 21)
        @test bounds(fccrop) == bounds(tablecrop) == ((-20, -5), (10, 30))
    end

end

using StableRNGs, StatsBase
test = rebuild(ga; name = :test)
@testset "sample" begin
    # test that all keywords work and return the same thing as extract
    @test all(Rasters.sample(StableRNG(123), test, 2) .=== extract(test, [(2.0,2.0), (1.0,2.0)]))
    @test all(Rasters.sample(StableRNG(123), st2, 2) .=== extract(st2, [(2,2), (1,2)]))
    @test all(Rasters.sample(StableRNG(123), test, 2; geometry = false) .=== extract(test, [(2.0,2.0), (1.0,2.0)]; geometry = false))
    @test all(Rasters.sample(StableRNG(123), st2, 2; geometry = false) .=== extract(st2, [(2,2), (1,2)]; geometry = false))
    
    @test all(Rasters.sample(StableRNG(123), test, 2, skipmissing = true) .=== extract(test, [(2.0,1.0), (2.0,1.0)], skipmissing = true))
    @test all(Rasters.sample(StableRNG(123), st2, 2, skipmissing = true) .=== extract(st2, [(2,1), (2,1)], skipmissing = true))

    @test all(
        Rasters.sample(StableRNG(123), test, 2, skipmissing = true, index = true) .=== 
        extract(test, [(2.0,1.0), (2.0,1.0)], skipmissing = true, index = true))
    @test all(
        Rasters.sample(StableRNG(123), st2, 2, skipmissing = true, index = true) .=== 
        extract(st2, [(2,1), (2,1)], skipmissing = true, index = true))

    kws = [(:geometry => false,), (:index => true,), ()]
    for kw in kws
        @test all(        
            Rasters.sample(StableRNG(123), test, 2; skipmissing = true, kw...) .=== 
            extract(test, [(2.0,1.0), (2.0,1.0)]; skipmissing = true, kw...)
            )
            
        @test all(
            Rasters.sample(StableRNG(123), st2, 2; skipmissing = true, kw...) .== 
            extract(st2, [(2,1), (2,1)]; skipmissing = true, kw...)
        )
    end
    @test all(Rasters.sample(StableRNG(123), st2, 2, name = (:a,)) .=== extract(st2, [(2,2), (1,2)], name = (:a,)))

    # extract and sample return the same type
    @test eltype(Rasters.sample(StableRNG(123), test, 2, index = true)) == eltype(extract(test, [(2.0,1.0), (2.0,1.0)], index = true))
    @test eltype(Rasters.sample(StableRNG(123), st2, 2, index = true)) == eltype(extract(st2, [(2,1), (2,1)], index = true))

    @test all(
        Rasters.sample(StableRNG(123), test, 2, weights = DimArray([1,1000], X(1:2)), skipmissing = true) .===
        [
            (geometry = (2.0,1.0), test = 2.0f0)
            (geometry = (2.0,1.0), test = 2.0f0)
        ]
    )

    @test all(
        Rasters.sample(StableRNG(123), test, 2, weights = DimArray([1,1000], X(1:2)), skipmissing = true, weightstype = StatsBase.FrequencyWeights) .===
        [
            (geometry = (2.0,1.0), test = 2.0f0)
            (geometry = (2.0,1.0), test = 2.0f0)
        ]
    )

    @test all(
        Rasters.sample(StableRNG(123), test, 2, skipmissing = true, replace = false) .===
        [
            (geometry = (2.0,1.0), test = 2.0f0)
            (geometry = (1.0,2.0), test = 7.0f0)
        ]
    )
    @test all(
        Rasters.sample(StableRNG(123), test, 2, skipmissing = true, replace = false, ordered = true) .===
        [
            (geometry = (2.0,1.0), test = 2.0f0)
            (geometry = (1.0,2.0), test = 7.0f0)
        ]
    )

    # test providing a geometry type works
    @test typeof(first(
        Rasters.sample(StableRNG(123), test, 2, geometry = (X = X, Y = Y))
    ).geometry) <: NamedTuple{(:X, :Y)}

    # test this works without providing n
    @test Rasters.sample(test, geometry = (X = X, Y = Y)) isa NamedTuple{(:geometry, :test)}

    @test_throws "strictly positive" Rasters.sample(StableRNG(123), test, 3, skipmissing = true, replace = false)
    @test_throws "Cannot draw" Rasters.sample(StableRNG(123), test, 5, replace = false)
end
