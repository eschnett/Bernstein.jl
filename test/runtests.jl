using Bernstein

using Random
using SimplexQuad
using StaticArrays
using Test

const Dmax = 5

const Rat128 = Rational{Int128}

# Random rationals
function Base.rand(rng::AbstractRNG, ::Random.SamplerType{Rational{T}}) where {T}
    return Rational{T}(T(rand(rng, -1000:1000)) // 1000)
end

# Check whether all matrix rows are unique
function unique_rows(A::AbstractArray{T,2}) where {T}
    nr = size(A, 1)
    for i in 1:nr, j in (i + 1):nr
        A[i, :] == A[j, :] && return false
    end
    return true
end

@testset "Barycentric coordinates for general simplices D=$D" for D in 0:Dmax
    T = Rat128
    N = D + 1

    for iter in 1:100
        s = rand(SMatrix{N,D,T})
        while !unique_rows(s)
            s = rand(SMatrix{N,D,T})
        end
        p = rand(SVector{D,T})
        λ = cartesian2barycentric(s, p)
        @test sum(λ) == 1
        setup = cartesian2barycentric_setup(s)
        λ′ = cartesian2barycentric(setup, p)
        @test λ′ == λ
        p′ = barycentric2cartesian(s, λ)
        @test p′ == p

        s2 = SVector{N,SVector{D,T}}(SVector{D,T}(s[i, a] for a in 1:D)
                                     for i in 1:N)
        p2 = p
        λ2 = cartesian2barycentric(s2, p2)
        @test λ2 == λ
        setup2 = cartesian2barycentric_setup(s2)
        λ2′ = cartesian2barycentric(setup2, p2)
        @test λ2′ == λ′
        p2′ = barycentric2cartesian(s2, λ2)
        @test p2′ == p′
    end
end

@testset "Simple Bernstein polynomials D=$D" for D in 0:Dmax
    T = Rat128
    N = D + 1

    s = SMatrix{N,D,T}(T(a + 1 == i) for i in 1:N, a in 1:D)

    if D == 0
        b0() = 1
        x = SVector{D,T}()
        @test bernstein(s, SVector(0), x) == b0()
    elseif D == 1
        b00(x) = oftype(x, 1)
        b01(x) = x
        b10(x) = 1 - x
        b02(x) = x^2
        b11(x) = 2 * x * (1 - x)
        b20(x) = (1 - x)^2
        b03(x) = x^3
        b12(x) = 3 * x^2 * (1 - x)
        b21(x) = 3 * x * (1 - x)^2
        b30(x) = (1 - x)^3
        for x1 in 0:(1 // 3):1
            x = SVector{D,T}(x1)
            @test bernstein(s, SVector(0, 0), x) == b00(x...)
            @test bernstein(s, SVector(0, 1), x) == b01(x...)
            @test bernstein(s, SVector(1, 0), x) == b10(x...)
            @test bernstein(s, SVector(0, 2), x) == b02(x...)
            @test bernstein(s, SVector(1, 1), x) == b11(x...)
            @test bernstein(s, SVector(2, 0), x) == b20(x...)
            @test bernstein(s, SVector(0, 3), x) == b03(x...)
            @test bernstein(s, SVector(1, 2), x) == b12(x...)
            @test bernstein(s, SVector(2, 1), x) == b21(x...)
            @test bernstein(s, SVector(3, 0), x) == b30(x...)
        end
    elseif D == 2
        b000(x, y) = oftype(x, 1)
        b001(x, y) = y
        b010(x, y) = x
        b100(x, y) = 1 - x - y
        b002(x, y) = y^2
        b011(x, y) = 2 * x * y
        b020(x, y) = x^2
        b101(x, y) = 2 * (1 - x - y) * y
        b110(x, y) = 2 * (1 - x - y) * x
        b200(x, y) = (1 - x - y)^2
        for x1 in 0:(1 // 4):(1 // 2), x2 in 0:(1 // 4):(1 // 2)
            x = SVector{D,T}(x1, x2)
            @test bernstein(s, SVector(0, 0, 0), x) == b000(x...)
            @test bernstein(s, SVector(0, 0, 1), x) == b001(x...)
            @test bernstein(s, SVector(0, 1, 0), x) == b010(x...)
            @test bernstein(s, SVector(1, 0, 0), x) == b100(x...)
            @test bernstein(s, SVector(0, 0, 2), x) == b002(x...)
            @test bernstein(s, SVector(0, 1, 1), x) == b011(x...)
            @test bernstein(s, SVector(0, 2, 0), x) == b020(x...)
            @test bernstein(s, SVector(1, 0, 1), x) == b101(x...)
            @test bernstein(s, SVector(1, 1, 0), x) == b110(x...)
            @test bernstein(s, SVector(2, 0, 0), x) == b200(x...)
        end
    elseif D == 3
        b0000(x, y, z) = oftype(x, 1)
        b0001(x, y, z) = z
        b0010(x, y, z) = y
        b0100(x, y, z) = x
        b1000(x, y, z) = 1 - x - y - z
        for x1 in 0:(1 // 2):(1 // 2),
x2 in 0:(1 // 2):(1 // 2),
x3 in 0:(1 // 2):(1 // 2)

            x = SVector{D,T}(x1, x2, x3)
            @test bernstein(s, SVector(0, 0, 0, 0), x) == b0000(x...)
            @test bernstein(s, SVector(0, 0, 0, 1), x) == b0001(x...)
            @test bernstein(s, SVector(0, 0, 1, 0), x) == b0010(x...)
            @test bernstein(s, SVector(0, 1, 0, 0), x) == b0100(x...)
            @test bernstein(s, SVector(1, 0, 0, 0), x) == b1000(x...)
        end
    end
end

# Overlap integrals between Bernstein polynomials (which are neither
# orthogonal nor normalized)
const bdbs = Dict{Tuple{Int,Int},Matrix{Float64}}()

@testset "Overlap between Bernstein polynomials D=$D P=$P" for D in 1:Dmax,
P in 1:min(2, 10 - 2D)

    T = Float64
    i2α, bdb = bernstein_products(T, D, P)
    bdbs[(D, P)] = bdb
    # println("|(D=$D,P=$P)|=$(length(bdb))")
    nα = binomial(P + D, D)
    @test size(bdb) == (nα, nα)
end

#TODO """
#TODO Integrate a function f over a simplex using Gauss quadrature
#TODO """
#TODO function integrate(f::F, X::SMatrix{D, N, T}, W::SVector{N, T}
#TODO                    ) where {F, D, N, T}
#TODO     sum(W[i] * f(X[:,i]) for i in 1:N)
#TODO end
#TODO 
#TODO """
#TODO Expand a function in P-th order Bernstein polynomials"
#TODO """
#TODO function expand(f, s::SMatrix{D, N, T}, P::Int) where {D, N, T}
#TODO     @assert N == D+1
#TODO 
#TODO     setup = bernstein_setup(s)
#TODO 
#TODO     P′ = P
#TODO     #TODO <https://github.com/JuliaArrays/StaticArrays.jl/issues/791>
#TODO     X1, W1 = simplexquad(P′, collect(s)')
#TODO     @assert size(X1) == (P^D, D)
#TODO     @assert size(W1) == (P^D,)
#TODO     X = SMatrix{D, P^D}(X1')
#TODO     W = SVector{P^D}(W1)
#TODO 
#TODO     αmin = CartesianIndex(ntuple(d -> 0, N))
#TODO     αmax = CartesianIndex(ntuple(d -> P, N))
#TODO     αs = SVector{N, Int}[]
#TODO     U = typeof(f(X[:,1]))
#TODO     cs = U[]
#TODO     for α1 in αmin:αmax
#TODO         α = SVector(α1.I)
#TODO         sum(α) == P || continue
#TODO 
#TODO         b(x) = bernstein(setup, α, x)
#TODO         k(x) = b(x)' * f(x)
#TODO         c = integrate(k, X, W)
#TODO 
#TODO         push!(αs, α)
#TODO         push!(cs, c)
#TODO     end
#TODO 
#TODO     # TODO: The Bernstein polynomials are not orthogonal!
#TODO     αs, cs
#TODO end
#TODO 
#TODO """
#TODO Evaluate a function given in Bernstein polynomial coefficients
#TODO """
#TODO function evaluate(s::SMatrix{D, N, T}, P::Int, cs::Vector{U}, x::SVector{D, T}
#TODO                   ) where {D, N, T, U}
#TODO     @assert N == D+1
#TODO 
#TODO     setup = bernstein_setup(s)
#TODO 
#TODO     αmin = CartesianIndex(ntuple(d -> 0, N))
#TODO     αmax = CartesianIndex(ntuple(d -> P, N))
#TODO     i = 0
#TODO     f = zero(U)
#TODO     for α1 in αmin:αmax
#TODO         α = SVector(α1.I)
#TODO         sum(α) == P || continue
#TODO 
#TODO         b(x) = bernstein(setup, α, x)
#TODO         i += 1
#TODO         c = cs[i]
#TODO 
#TODO         f += c * b(x)
#TODO     end
#TODO     @assert i == length(cs)
#TODO 
#TODO     f
#TODO end
#TODO 
#TODO 
#TODO 
#TODO #TODO @testset "Expand and evaluate some polynomials D=$D P=$P" for D in 1:Dmax, P in 1:11-2D
#TODO @testset "Expand and evaluate some polynomials D=$D P=$P" for D in 1:1, P in 1:1
#TODO     T = Float64
#TODO     N = D+1
#TODO 
#TODO     s = SMatrix{D, N}(T(a+1 == i) for a in 1:D, i in 1:N)
#TODO 
#TODO     f(x) = sum(x[i] for i in 1:D)
#TODO     αs, cs = expand(f, s, P)
#TODO     @show αs cs
#TODO 
#TODO     @error "take non-orthonormality into account"
#TODO 
#TODO     if D == 1
#TODO         for x1 in 0:1//P:1
#TODO             x = SVector{D, T}(x1)
#TODO             fx = evaluate(s, P, cs, x)
#TODO             @show x f(x) fx fx-f(x)
#TODO         end
#TODO     end
#TODO end
