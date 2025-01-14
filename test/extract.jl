using Rasters, Test, DataFrames, Extents
import GeoInterface as GI
using Rasters.Lookups, Rasters.Dimensions 

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

dimz = (X(9.0:1.0:10.0), Y(0.1:0.1:0.2))
rast = Raster(Union{Int,Missing}[1 2; 3 4], dimz; name=:test, missingval=missing)
rast2 = Raster([5 6; 7 8], dimz; name=:test2, missingval=5)
rast_m = Raster([1 2; 3 missing], dimz; name=:test, missingval=missing)
st = RasterStack(rast, rast2)

pts = [missing, (9.0, 0.1), (9.0, 0.2), (10.0, 0.3), (10.0, 0.2)]
poly = GI.Polygon([[(8.0, 0.0), (11.0, 0.0), (11.0, 0.4), (8.0, 0.0)]])
linestring = GI.LineString([(8.0, 0.0), (9.5, 0.0), (10.0, 0.4)])
line = GI.Line([(8.0, 0.0), (12.0, 0.4)])
table = (geometry=pts, foo=zeros(4))

@testset "Points" begin
    @testset "From Raster" begin
        @testset "skipmissing=false" begin
            ex = extract(rast, pts; skipmissing=false)
            T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},test::Union{Missing,Int64}}
            @test eltype(ex) == T
            @test all(ex .=== T[
                (geometry = missing, test = missing)
                (geometry = (9.0, 0.1), test=1)
                (geometry = (9.0, 0.2), test=2)
                (geometry = (10.0, 0.3), test=missing)
                (geometry = (10.0, 0.2), test=4)
            ])
        end
        @testset "skipmissing=true" begin
            ex = extract(rast_m, pts; skipmissing=true);
            T = @NamedTuple{geometry::Tuple{Float64, Float64}, test::Int64}
            @test eltype(ex) == T
            @test all(ex .=== T[(geometry = (9.0, 0.1), test = 1), (geometry = (9.0, 0.2), test = 2)])
        end

        @testset "skipmissing=true, geometry=false" begin
            ex = extract(rast_m, pts; skipmissing=true, geometry=false)
            T = @NamedTuple{test::Int64}
            @test eltype(ex) == T
            @test all(ex .=== T[(test = 1,), (test = 2,)])
            @test all(extract(rast_m, pts; skipmissing=true, geometry=false, index=true) .=== [
                (index = (1, 1), test = 1,)
                (index = (1, 2), test = 2,)
            ])
        end
        @testset "reverse points" begin
            # NamedTuple (reversed) points - tests a Table that iterates over points
            T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},test::Union{Missing,Int64}}
            @test all(extract(rast, [(Y=0.1, X=9.0), (Y=0.2, X=10.0), (Y=0.3, X=10.0)]) .=== T[
                (geometry = (9.0,  0.1), test = 1)
                (geometry = (10.0, 0.2), test = 4)
                (geometry = (10.0, 0.3), test = missing)
            ])
        end
        @testset "Vector points" begin
            @test extract(rast, [[9.0, 0.1], [10.0, 0.2]]) == [
                (geometry = (9.0, 0.1), test = 1)
                (geometry = (10.0, 0.2), test = 4)
            ]
        end
    end

    @testset "From RasterStack" begin
        T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},test::Union{Missing,Int64},test2::Union{Missing,Int64}}
        @test all(extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]) .=== T[
            (geometry = missing, test = missing, test2 = missing)
            (geometry = (9.0, 0.1), test = 1, test2 = 5)
            (geometry = (10.0, 0.2), test = 4, test2 = 8)
            (geometry = (10.0, 0.3), test = missing, test2 = missing)
        ])
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; skipmissing=true) == [
            (geometry = (10.0, 0.2), test = 4, test2 = 8)
        ]
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; skipmissing=true, geometry=false) == [
            (test = 4, test2 = 8)
        ]
        T = @NamedTuple{index::Union{Missing, Tuple{Int,Int}}, test::Union{Missing, Int64}, test2::Union{Missing, Int64}}
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; skipmissing=true, geometry=false, index=true) == T[
            (index = (2, 2), test = 4, test2 = 8)
        ]
        # Subset with `names`
        T = @NamedTuple{geometry::Union{Missing, Tuple{Float64, Float64}}, test2::Union{Missing, Int64}}
        @test all(extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; name=(:test2,)) .=== T[
            (geometry = missing, test2 = missing)
            (geometry = (9.0, 0.1), test2 = 5)
            (geometry = (10.0, 0.2), test2 = 8)
            (geometry = (10.0, 0.3), test2 = missing)
        ])
        # Subset with `names` and `skipmissing` with mixed missingvals
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; name=(:test2,), skipmissing=true) == [
            (geometry = (10.0, 0.2), test2 = 8)
        ]
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; name=(:test,), skipmissing=true) == [
            (geometry = (9.0, 0.1), test = 1)
            (geometry = (10.0, 0.2), test = 4)
        ]
    end

    @testset "Errors" begin
        # Missing coord errors
        @test_throws ArgumentError extract(rast, [(0.0, missing), (9.0, 0.1), (9.0, 0.2), (10.0, 0.3)])
        @test_throws ArgumentError extract(rast, [(9.0, 0.1), (0.0, missing), (9.0, 0.2), (10.0, 0.3)])
        @test_throws ArgumentError extract(rast, [(X=0.0, Y=missing), (9.0, 0.1), (9.0, 0.2), (10.0, 0.3)])
    end
end

@testset "Polygons" begin
    # Extract a polygon
    T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
    ])
    # Test all the keyword combinations
    @test all(extract(rast_m, poly) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
    ])
    T = @NamedTuple{test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly; geometry=false) .=== T[
        (test = 1,)
        (test = 3,)
        (test = missing,)
    ])
    T = @NamedTuple{index::Union{Missing,Tuple{Int,Int}},test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly; geometry=false, index=true) .=== T[
        (index = (1, 1), test = 1)
        (index = (2, 1), test = 3)
        (index = (2, 2), test = missing)
    ])
    T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},index::Union{Missing,Tuple{Int,Int}},test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly; index=true) .=== T[
         (geometry = (9.0, 0.1), index = (1, 1), test = 1)
         (geometry = (10.0, 0.1), index = (2, 1), test = 3)
         (geometry = (10.0, 0.2), index = (2, 2), test = missing)
    ])
    @test extract(rast_m, poly; skipmissing=true) == [
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
    ]                                                         
    @test extract(rast_m, poly; skipmissing=true, geometry=false) == [
        (test = 1,)
        (test = 3,)
    ]                                                         
    @test extract(rast_m, poly; skipmissing=true, geometry=false, index=true) == [
        (index = (1, 1), test = 1)
        (index = (2, 1), test = 3)
    ]                                                         
    @test extract(rast_m, poly; skipmissing=true, index=true) == [
        (geometry = (9.0, 0.1), index = (1, 1), test = 1)
        (geometry = (10.0, 0.1), index = (2, 1), test = 3)
    ]          
    @test extract(rast2, poly; skipmissing=true) == [
        (geometry = (10.0, 0.1), test2 = 7)
        (geometry = (10.0, 0.2), test2 = 8)
    ]                                               
    T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
    ])
    @test extract(rast_m, poly; skipmissing=true, geometry=false, id=true) == [
        (id=1, test = 1,)
        (id=1, test = 3,)
    ]                                                         
    @test extract(rast_m, [poly, poly]; skipmissing=true, geometry=false, id=true) == [
        (id=1, test = 1,)
        (id=1, test = 3,)
        (id=2, test = 1,)
        (id=2, test = 3,)
    ]                                                         
    @test extract(rast_m, [poly, poly]; skipmissing=true, geometry=false, id=true, flatten=false) == [
        [(id=1, test = 1,), (id=1, test = 3,)],
        [(id=2, test = 1,), (id=2, test = 3,)],
    ]                                                         

    @testset "Vector of polygons" begin
        ex = extract(rast_m, [poly, poly])
        @test eltype(ex) == T
        @test all(ex .=== T[
            (geometry = (9.0, 0.1), test = 1)
            (geometry = (10.0, 0.1), test = 3)
            (geometry = (10.0, 0.2), test = missing)
            (geometry = (9.0, 0.1), test = 1)
            (geometry = (10.0, 0.1), test = 3)
            (geometry = (10.0, 0.2), test = missing)
        ])
    end
end

@testset "Extract a linestring" begin
    T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},test::Union{Missing,Int64}}
    Tsm = @NamedTuple{geometry::Tuple{Float64,Float64},test::Int64}
    linestrings = [linestring, linestring, linestring]
    fc = GI.FeatureCollection(map(GI.Feature, linestrings))

    # Single LineString
    @test extract(rast, linestring) isa Vector{T}
    @test all(extract(rast_m, linestring) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
    ])
    @test all(extract(rast_m, linestring; skipmissing=true) .=== Tsm[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
    ])

    # Multiple linstrings
    # Test all combinations of skipmissing, flatten and threaded
    @test extract(rast_m, linestrings; skipmissing=true) isa Vector{Tsm}
    @test extract(rast_m, fc; skipmissing=true) == 
          extract(rast_m, fc; skipmissing=true, threaded=true) == 
          extract(rast_m, linestrings; skipmissing=true) == 
          extract(rast_m, linestrings; skipmissing=true, threaded=true) == Tsm[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
    ]

    @test extract(rast_m, linestrings; skipmissing=false) isa Vector{T}
    @test all(
              extract(rast_m, fc; skipmissing=false) .=== 
              extract(rast_m, fc; skipmissing=false, threaded=true) .===
              extract(rast_m, linestrings; skipmissing=false) .=== 
              extract(rast_m, linestrings; skipmissing=false, threaded=true) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
    ])

    @test extract(rast_m, linestrings; skipmissing=true, flatten=false) isa Vector{Vector{Tsm}}
    @test extract(rast_m, linestrings; skipmissing=true, flatten=false) == 
        extract(rast_m, linestrings; skipmissing=true, flatten=false, threaded=true) == Vector{Tsm}[
        [(geometry = (9.0, 0.1), test = 1), (geometry = (10.0, 0.1), test = 3)],
        [(geometry = (9.0, 0.1), test = 1), (geometry = (10.0, 0.1), test = 3)],
        [(geometry = (9.0, 0.1), test = 1), (geometry = (10.0, 0.1), test = 3)],
    ]

    # Nested Vector holding missing needs special handling to check equality
    @test extract(rast_m, linestrings; skipmissing=false, flatten=false) isa Vector{Vector{T}}
    ref = Vector{T}[
        [(geometry = (9.0, 0.1), test = 1), (geometry = (10.0, 0.1), test = 3), (geometry = (10.0, 0.2), test = missing)] ,
        [(geometry = (9.0, 0.1), test = 1), (geometry = (10.0, 0.1), test = 3), (geometry = (10.0, 0.2), test = missing)],
        [(geometry = (9.0, 0.1), test = 1), (geometry = (10.0, 0.1), test = 3), (geometry = (10.0, 0.2), test = missing)],
    ]
    matching(a, b) = all(map(===, a, b))
    @test all(map(matching, extract(rast_m, linestrings; skipmissing=false, flatten=false, threaded=false), ref))
    @test all(map(matching, extract(rast_m, linestrings; skipmissing=false, flatten=false, threaded=true), ref))
end

@testset "Extract a line" begin
    T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},test::Union{Missing,Int64}}
    Tsm = @NamedTuple{geometry::Tuple{Float64,Float64},test::Int64}
    lines = [line, line]
    fc = GI.FeatureCollection(map(GI.Feature, lines))

    # Single LineString
    @test extract(rast, line) isa Vector{T}
    @test all(extract(rast_m, line) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
        (geometry = (10.0, 0.2), test = missing)
    ])
    @test all(extract(rast_m, line; skipmissing=true) .=== Tsm[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
    ])

    # Multiple linstrings
    # Test all combinations of skipmissing, flatten and threaded
    @test extract(rast_m, lines; skipmissing=true) isa Vector{Tsm}
    @test extract(rast_m, fc; skipmissing=true) == 
          extract(rast_m, fc; skipmissing=true, threaded=true) == 
          extract(rast_m, lines; skipmissing=true) == 
          extract(rast_m, lines; skipmissing=true, threaded=true) == Tsm[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
    ]

    @test extract(rast_m, lines; skipmissing=false) isa Vector{T}
    @test all(
              extract(rast_m, fc; skipmissing=false) .=== 
              extract(rast_m, fc; skipmissing=false, threaded=true) .===
              extract(rast_m, lines; skipmissing=false) .=== 
              extract(rast_m, lines; skipmissing=false, threaded=true) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
        (geometry = (10.0, 0.2), test = missing)
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
        (geometry = (10.0, 0.2), test = missing)
    ])

    Tsm_i = @NamedTuple{id::Int,geometry::Tuple{Float64,Float64},test::Int64}
    @test extract(rast, lines; skipmissing=true, id=true) isa Vector{Tsm_i}
    @test all(extract(rast_m, fc; skipmissing=true, id=true) .=== 
              extract(rast_m, fc; skipmissing=true, threaded=true, id=true) .===
              extract(rast_m, lines; skipmissing=true, id=true) .=== 
              extract(rast_m, lines; skipmissing=true, threaded=true, id=true) .=== Tsm_i[
        (id=1, geometry = (9.0, 0.1), test = 1)
        (id=1, geometry = (9.0, 0.2), test = 2)
        (id=2, geometry = (9.0, 0.1), test = 1)
        (id=2, geometry = (9.0, 0.2), test = 2)
    ])

    @test extract(rast_m, lines; skipmissing=true, flatten=false) isa Vector{Vector{Tsm}}
    @test extract(rast_m, lines; skipmissing=true, flatten=false) == 
        extract(rast_m, lines; skipmissing=true, flatten=false, threaded=true) == Vector{Tsm}[
        [(geometry = (9.0, 0.1), test = 1), (geometry = (9.0, 0.2), test = 2)],
        [(geometry = (9.0, 0.1), test = 1), (geometry = (9.0, 0.2), test = 2)],
    ]

    # Nested Vector holding missing needs special handling to check equality
    @test extract(rast_m, lines; skipmissing=false, flatten=false) isa Vector{Vector{T}}
    ref = Vector{T}[
        [(geometry = (9.0, 0.1), test = 1), (geometry = (9.0, 0.2), test = 2), (geometry = (10.0, 0.2), test = missing)] ,
        [(geometry = (9.0, 0.1), test = 1), (geometry = (9.0, 0.2), test = 2), (geometry = (10.0, 0.2), test = missing)],
    ]
    matching(a, b) = all(map(===, a, b))
    @test all(map(matching, extract(rast_m, lines; skipmissing=false, flatten=false, threaded=false), ref))
    @test all(map(matching, extract(rast_m, lines; skipmissing=false, flatten=false, threaded=true), ref))
end

@testset "Table" begin
    T = @NamedTuple{geometry::Union{Missing,Tuple{Float64, Float64}}, test::Union{Missing, Int64}}
    @test all(extract(rast, table) .=== T[
        (geometry = missing, test = missing)
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
        (geometry = (10.0, 0.3), test = missing)
        (geometry = (10.0, 0.2), test = 4)
    ])
    @test extract(rast, table; skipmissing=true) == [
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
        (geometry = (10.0, 0.2), test = 4)
    ]
    @test extract(rast, table; skipmissing=true, geometry=false) == [
        (test = 1,)
        (test = 2,)
        (test = 4,)
    ]
    @test extract(rast, table; skipmissing=true, geometry=false, index=true) == [
        (index = (1, 1), test = 1,)
        (index = (1, 2), test = 2,)
        (index = (2, 2), test = 4,)
    ]
    @test extract(rast, table; skipmissing=true, geometry=false, id=true) == [
        (id=2, test = 1,)
        (id=3, test = 2,)
        (id=5, test = 4,)
    ]
    T = @NamedTuple{id::Int, test::Union{Missing, Int64}}
    @test all(extract(rast, table; skipmissing=false, geometry=false, id=true) .=== T[
        (id=1, test = missing,)
        (id=2, test = 1,)
        (id=3, test = 2,)
        (id=4, test = missing,)
        (id=5, test = 4,)
    ])
    @test_throws ArgumentError extract(rast, (foo = zeros(4),))
end

@testset "Empty geoms" begin
    @test all(extract(rast, []) .=== @NamedTuple{geometry::Missing,test::Union{Missing,Int}}[])
    @test all(extract(rast, []; geometry=false) .=== @NamedTuple{test::Tuple{Missing}}[])
    @test all(extract(rast, [missing, missing]; geometry=false) .=== @NamedTuple{test::Union{Missing,Int}}[
         (test = missing,)
         (test = missing,)
    ])
    @test typeof(extract(rast, []; geometry=false, skipmissing=true)) == Vector{@NamedTuple{test::Int}}
    @test typeof(extract(rast, [missing, missing]; geometry=false, skipmissing=true)) == Vector{@NamedTuple{test::Int}}
end

#=
Some benchmark and plotting code

using BenchmarkTools
using ProfileView
using Chairmarks
using Cthulhu
using Rasters
using DataFrames
import GeoInterface as GI
using Rasters.Lookups, Rasters.Dimensions 

dimz = (X(5.0:0.1:15.0; sampling=Intervals(Start())), Y(-0.1:0.01:0.5; sampling=Intervals(Start())))
rast_m = Raster(rand([missing, 1, 2], dimz); name=:test, missingval=missing)
rast = Raster(rand(Int, dimz); name=:test2, missingval=5)
st = RasterStack((a=rast, b=rast, c=Float32.(rast)))
ks = ntuple(x -> gensym(), 100)
st100 = RasterStack(NamedTuple{ks}(ntuple(x -> rast, length(ks))))

poly = GI.Polygon([[(8.0, 0.0), (11.0, 0.0), (11.0, 0.4), (8.0, 0.0)]])
linestring = GI.LineString([(8.0, 0.0), (11.0, 0.0), (11.0, 0.4), (8.0, 0.0)])
line = GI.Line([(8.0, 0.0), (12.0, 4.0)])
polys = [poly for i in 1:10000]
lines = [line for i in 1:100000]
linestrings = [linestring for i in 1:100000]

extract(rast, linestring)
extract(rast, polys; threaded=false)
@time extract(rast, polys; geometry=false, threaded=true, flatten=true) |> length
@time extract(rast, lines; geometry=false, threaded=true, flatten=true) |> length
@time extract(st, polys; geometry=false, threaded=true, flatten=true) |> length
@time extract(st, lines; geometry=false, threaded=false, flatten=true) |> length

@time extract(rast_m, lines; geometry=false, threaded=true, flatten=false) |> length
@time extract(rast, lines; geometry=false, threaded=true, flatten=false) |> length
@time extract(rast_m, lines; skipmissing=true, geometry=false, threaded=true, flatten=false) |> length
@time extract(rast,   lines; skipmissing=true, geometry=false, threaded=true, flatten=false) |> length
@time extract(st,     lines; skipmissing=true, geometry=false, threaded=true, flatten=false) |> length
@time extract(st100, lines; geometry=false, threaded=false, flatten=false) |> length
@time extract(st100, lines; geometry=false, threaded=true, flatten=false) |> length
@time extract(rast, lines; geometry=false, threaded=true, flatten=true) |> length
@time extract(rast, lines; threaded=true, flatten=true, progress=false) |> length
@time extract(rast, lines; threaded=false, flatten=false, progress=false) |> length
@benchmark extract(rast, lines; threaded=false, flatten=true)
@benchmark extract(rast, lines; threaded=true, flatten=true)
@benchmark extract(rast, lines; threaded=false, flatten=false)
@benchmark extract(rast, lines; threaded=true, flatten=false)
@time extract(rast, polys; flatten=false);
extract(rast, polys; flatten=true);
extract(st100, linestring)
extract(st100, poly)

@profview extract(rast, lines)
@profview extract(rast, lines; threaded=true)
@profview extract(rast, linestrings)
@profview extract(rast_m, linestrings)
@profview extract(rast, polys)
@profview extract(rast, polys; flatten=false)
@profview extract(st, linestring)
@profview extract(st, poly)

# Profile running one function a lot. 
# People will do this
f(rast, ls, n; skipmissing=true) = for _ in 1:n extract(rast, ls; skipmissing) end
@profview f(rast, polies, 10; skipmissing=false)
@profview f(rast, poly, 10000; skipmissing=false)
@profview f(rast_m, poly, 10000; skipmissing=false)
@profview f(rast, linestring, 10000; skipmissing=false)
@profview f(rast_m, linestring, 10000; skipmissing=false)
@profview f(st, linestring, 10000; skipmissing=false)
@profview f(rast, linestring, 10000)
@profview f(st, linestring, 100000; skipmissing=false)
@profview f(st, linestring, 100000; skipmissing=true)
@profview f(st100, linestring, 1000)
@profview f(st100, line, 1000)
@profview f(st100, poly, 10000; skipmissing=true)

@profview f(st100, linestring, 1000; skipmissing=true)
@profview f(st100, linestring, 1000; skipmissing=false)
@profview f(st100, linestring, 1000; names=skipmissing=false)
keys(st100)[1:50]



# Demo plots
using GLMakie
using GeoInterfaceMakie

Makie.plot(rast; colormap=:Reds)
Makie.plot!(rast[Touches(GI.extent(first(polys)))] ; colormap=:Greens)
Makie.plot!(polys; alpha=0.2)
ps = getindex.(extract(rast, poly; index=true, flatten=true), :geometry)
Makie.scatter!(ps; color=:yellow)

Makie.plot(rast; colormap=:Reds)
Makie.plot!(rast[Touches(GI.extent(first(polys)))] ; colormap=:Greens)
Makie.plot!(poly; alpha=0.7)
ps = getindex.(extract(rast, poly; index=true, flatten=true, boundary=:touches), :geometry);
Makie.scatter!(ps; color=:pink)

Makie.plot(rast; colormap=:Reds)
Makie.plot!(rast[Touches(GI.extent(first(polys)))] ; colormap=:Greens)
Makie.plot!(polys; alpha=0.7)
ps = getindex.(extract(rast, poly; index=true, flatten=true, boundary=:inside), :geometry)
Makie.scatter!(ps; color=:pink)

Makie.plot(rast; colormap=:Reds)
Makie.plot!(rast[Touches(GI.extent(first(polys)))] ; colormap=:Greens)
Makie.plot!(linestring; color=:violet, linewidth=5)
ps = getindex.(extract(rast, linestring; index=true), :geometry);
Makie.scatter!(ps; color=:yellow)

=#
