using Bernstein

using Grassmann
using Random
using StaticArrays
using Test



const Dmax = 5

const Rat128 = Rational{Int128}

# Random rationals
Base.rand(rng::AbstractRNG, ::Random.SamplerType{Rational{T}}) where {T} =
    Rational{T}(T(rand(rng, -1000:1000)) // 1000)



@testset "Barycentric coordinates for unit simplices D=$D" for D in 1:Dmax
    T = Rat128
    S = Signature(D)
    V = SubManifold(S)

    # Test simple algorithm
    for i in 1:100
        p = rand(Chain{V, 1, T})
        λ = cartesian2barycentric(p)
        @test sum(λ.v) == 1
        p′ = barycentric2cartesian(λ)
        @test p == p′
    end

    # Test generic algorithm
    N = D+1
    #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
    s = SVector{N}([Chain{V, 1}(SVector{D}([T(i+1==n) for i in 1:D]))
                    for n in 1:N])
    for i in 1:100
        p = rand(Chain{V, 1, T})
        λ = cartesian2barycentric(p)
        λ′ = cartesian2barycentric(s, p)
        @test λ == λ′
        p′ = barycentric2cartesian(s, λ)
        @test p == p′
    end
end

@testset "Barycentric coordinates for general simplices D=$D" for D in 1:Dmax
    T = Rat128
    S = Signature(D)
    V = SubManifold(S)
    N = D+1

    for i in 1:100
        #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
        s = SVector{N}([rand(Chain{V, 1, T}) for n in 1:N])
        p = rand(Chain{V, 1, T})
        λ = cartesian2barycentric(p)
        @test sum(λ.v) == 1
        p′ = barycentric2cartesian(λ)
        @test p == p′
    end
end



@testset "Simple Bernstein polynomials D=$D" for D in 1:Dmax
    T = Rat128
    S = Signature(D)
    V = SubManifold(S)
    N = D+1

    #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
    s = SVector{N}([Chain{V, 1}(SVector{D}([T(i+1==n) for i in 1:D]))
                    for n in 1:N])

    if D == 1
        b00(x) = oftype(x, 1)
        b01(x) = x
        b10(x) = 1-x
        b02(x) = x^2
        b11(x) = 2*x*(1-x)
        b20(x) = (1-x)^2
        b03(x) = x^3
        b12(x) = 3*x^2*(1-x)
        b21(x) = 3*x*(1-x)^2
        b30(x) = (1-x)^3
        for x1 in 0:1//3:1
            x = Chain{V, 1}(T(x1))
            @test bernstein(s, SVector(0,0), x) == b00(x.v...)
            @test bernstein(s, SVector(0,1), x) == b01(x.v...)
            @test bernstein(s, SVector(1,0), x) == b10(x.v...)
            @test bernstein(s, SVector(0,2), x) == b02(x.v...)
            @test bernstein(s, SVector(1,1), x) == b11(x.v...)
            @test bernstein(s, SVector(2,0), x) == b20(x.v...)
            @test bernstein(s, SVector(0,3), x) == b03(x.v...)
            @test bernstein(s, SVector(1,2), x) == b12(x.v...)
            @test bernstein(s, SVector(2,1), x) == b21(x.v...)
            @test bernstein(s, SVector(3,0), x) == b30(x.v...)
        end
    elseif D == 2
        b000(x,y) = oftype(x, 1)
        b001(x,y) = y
        b010(x,y) = x 
        b100(x,y) = 1 - x - y
        b002(x,y) = y^2
        b011(x,y) = 2 * x * y
        b020(x,y) = x^2
        b101(x,y) = 2 * (1 - x - y) * y
        b110(x,y) = 2 * (1 - x - y) * x
        b200(x,y) = (1 - x - y)^2
        for x1 in 0:1//4:1//2, x2 in 0:1//4:1//2
            x = Chain{V, 1}(T(x1), T(x2))
            @test bernstein(s, SVector(0,0,0), x) == b000(x.v...)
            @test bernstein(s, SVector(0,0,1), x) == b001(x.v...)
            @test bernstein(s, SVector(0,1,0), x) == b010(x.v...)
            @test bernstein(s, SVector(1,0,0), x) == b100(x.v...)
            @test bernstein(s, SVector(0,0,2), x) == b002(x.v...)
            @test bernstein(s, SVector(0,1,1), x) == b011(x.v...)
            @test bernstein(s, SVector(0,2,0), x) == b020(x.v...)
            @test bernstein(s, SVector(1,0,1), x) == b101(x.v...)
            @test bernstein(s, SVector(1,1,0), x) == b110(x.v...)
            @test bernstein(s, SVector(2,0,0), x) == b200(x.v...)
       end
    elseif D == 3
        b0000(x,y,z) = oftype(x, 1)
        b0001(x,y,z) = z
        b0010(x,y,z) = y
        b0100(x,y,z) = x 
        b1000(x,y,z) = 1 - x - y - z
        for x1 in 0:1//2:1//2, x2 in 0:1//2:1//2, x3 in 0:1//2:1//2
            x = Chain{V, 1}(T(x1), T(x2), T(x3))
            @test bernstein(s, SVector(0,0,0,0), x) == b0000(x.v...)
            @test bernstein(s, SVector(0,0,0,1), x) == b0001(x.v...)
            @test bernstein(s, SVector(0,0,1,0), x) == b0010(x.v...)
            @test bernstein(s, SVector(0,1,0,0), x) == b0100(x.v...)
            @test bernstein(s, SVector(1,0,0,0), x) == b1000(x.v...)
        end
    end
end
