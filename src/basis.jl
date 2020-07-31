# Bernstein basis

export bernstein_products
"""
Overlap integrals between Bernstein polynomials (which are neither
orthogonal nor normalized)
"""
function bernstein_products(::Type{T}, ::Val{D}, ::Val{P}) where {T,D,P}
    N = D + 1

    s = SMatrix{N,D}(T(a + 1 == i) for i in 1:N, a in 1:D)

    NP = P                # number of integration points per dimension
    NC = NP^D             # total number of integration points
    X1, W1 = simplexquad(NP, collect(s))
    @assert size(X1) == (NC, D)
    @assert size(W1) == (NC,)
    X = [SVector{D,T}(X1[n, a] for a in 1:D) for n in 1:NC]
    X::Vector{SVector{D,T}}
    # dump(X)
    W = W1

    setup = bernstein_setup(s)
    Λ = [cartesian2barycentric(setup, x) for x in X]
    Λ::Vector{SVector{N,T}}

    nα = binomial(P + D, D)
    bdb = Array{T}(undef, nα, nα)
    i2α = Array{SVector{N,Int}}(undef, nα)

    # Instead of re-evaluting the Bernstein polynomials at each
    # integration point, we could pre-calculate them. This would
    # increase storage significantly, but reduce computation time
    # significantly.

    αmin = CartesianIndex(ntuple(d -> 0, N))
    αmax = CartesianIndex(ntuple(d -> P, N))

    i = 0
    for αiI in αmin:αmax
        αi = SVector(αiI.I)
        sum(αi) == P || continue
        i += 1

        i2α[i] = αi

        j = 0
        for αjI in αmin:αmax
            αj = SVector(αjI.I)
            sum(αj) == P || continue
            j += 1

            αj >= αi || continue

            f(λ) = bernstein(αi, λ) * bernstein(αj, λ)
            bdb[i, j] = integrate_λ(f, Λ, W)
        end
        @assert j == nα
    end
    @assert i == nα

    return i2α, Symmetric(bdb)
end

bernstein_products(T, D, P) = bernstein_products(T, Val(D), Val(P))

@fastmath function integrate_λ(f::F, Λ::Vector{SVector{N,T}},
                               W::Vector{T}) where {F,N,T}
    @assert length(Λ) == length(W) > 0
    return @inbounds begin
        s = zero(W[1] * f(Λ[1]))
        for (λ, w) in Iterators.zip(Λ, W)
            s += w * f(λ)
        end
        return s
    end
end
