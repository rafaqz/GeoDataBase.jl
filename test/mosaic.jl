using Rasters, Dates, Statistics, Test
using Rasters.Lookups, Rasters.Dimensions 

@testset "mosaic" begin
    reg1 = Raster([0.1 0.2; 0.3 0.4], (X(2.0:1.0:3.0), Y(5.0:1.0:6.0)))
    reg2 = Raster([1.1 1.2; 1.3 1.4], (X(3.0:1.0:4.0), Y(6.0:1.0:7.0)))
    irreg1 = Raster([0.1 0.2; 0.3 0.4], (X([2.0, 3.0]), Y([5.0, 6.0])))
    irreg2 = Raster([1.1 1.2; 1.3 1.4], (X([3.0, 4.0]), Y([6.0, 7.0])))

    span_x1 = Explicit(vcat((1.5:1.0:2.5)', (2.5:1.0:3.5)'))
    span_x2 = Explicit(vcat((2.5:1.0:3.5)', (3.5:1.0:4.5)'))
    exp1 = Raster([0.1 0.2; 0.3 0.4], (X(Sampled([2.0, 3.0]; span=span_x1)), Y([5.0, 6.0])))
    exp2 = Raster([1.1 1.2; 1.3 1.4], (X(Sampled([3.0, 4.0]; span=span_x2)), Y([6.0, 7.0])))
    @test val(span(mosaic(first, exp1, exp2), X)) == [1.5 2.5 3.5; 2.5 3.5 4.5]
    @test all(mosaic(first, [reg1, reg2]) .=== 
              mosaic(first, irreg1, irreg2) .===
              mosaic(first, (irreg1, irreg2)) .=== 
              [0.1 0.2 missing; 
               0.3 0.4 1.2; 
               missing 1.3 1.4])
    @test all(mosaic(last, reg1, reg2) .===
              mosaic(last, irreg1, irreg2) .===
              mosaic(last, exp1, exp2) .=== [0.1 0.2 missing; 
                                             0.3 1.1 1.2; 
                                             missing 1.3 1.4])

    @test all(mosaic(first, [reverse(reg2; dims=Y), reverse(reg1; dims=Y)]) .=== 
        [missing 0.2 0.1; 
         1.2 1.1 0.3; 
         1.4 1.3 missing]
    )
            
    @testset "Generic functions" begin
        @test all(mosaic(xs -> count(x -> x > 0, xs), reg1, reg2) .===
            [1 1 missing
             1 2 1 
             missing 1 1]
        )
    end

   @testset "3 dimensions" begin
        A1 = Raster(ones(2, 2, 2), (X(2.0:-1.0:1.0), Y(5.0:1.0:6.0), Ti(DateTime(2001):Year(1):DateTime(2002))))
        A2 = Raster(zeros(2, 2, 2), (X(3.0:-1.0:2.0), Y(4.0:1.0:5.0), Ti(DateTime(2002):Year(1):DateTime(2003))))
        mean_mos = cat([missing missing missing
                    missing 1.0     1.0
                    missing 1.0     1.0    ],
                    [0.0     0.0 missing
                    0.0     0.5     1.0   
                    missing 1.0     1.0    ],
                    [0.0     0.0     missing
                    0.0     0.0     missing    
                    missing missing missing], dims=3)
        @test all(mosaic(mean, A1, A2) .=== 
                mosaic(mean, RasterStack(A1), RasterStack(A2)).layer1 .===
                mean_mos)
        @test mosaic(length, A1, A2; missingval=0) == cat([0 0 0
            0 1 1
            0 1 1],
            [1 1 0
            1 2 1
            0 1 1], 
            [1 1 0
            1 1 0
            0 0 0], dims=3)
    end
end
