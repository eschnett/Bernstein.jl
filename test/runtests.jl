using Bernstein

using ComputedFieldTypes
using Grassmann
using StaticArrays
using Test



const Dmax = 5

@testset "Barycentric coordinates for unit simplices D=$D" for D in 1:Dmax
    T = Float64
    S = Signature(D)
    V = SubManifold(S)

    # Test simple algorithm
    for i in 1:100
        p = rand(Chain{V, 1, T})
        λ = cartesian2barycentric(p)
        @test sum(λ.v) ≈ 1
        p′ = barycentric2cartesian(λ)
        @test p ≈ p′
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
        @test λ ≈ λ′
        p′ = barycentric2cartesian(s, λ)
        @test p ≈ p′
    end
end

@testset "Barycentric coordinates for general simplices D=$D" for D in 1:Dmax
    T = Float64
    S = Signature(D)
    V = SubManifold(S)
    N = D+1

    for i in 1:100
        #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
        s = SVector{N}([rand(Chain{V, 1, T}) for n in 1:N])
        p = rand(Chain{V, 1, T})
        λ = cartesian2barycentric(p)
        @test sum(λ.v) ≈ 1
        p′ = barycentric2cartesian(λ)
        @test p ≈ p′
    end
end
