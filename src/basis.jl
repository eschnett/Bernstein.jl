# Bernstein basis

export bernstein_products
"""
Overlap integrals between Bernstein polynomials (which are neither
orthogonal nor normalized)
"""
function bernstein_products(::Type{T}, ::Val{D}, ::Val{P}) where {T, D, P}
    N = D+1

    s = SMatrix{D, N}(T(a+1 == i) for a in 1:D, i in 1:N)

    NP = P                # number of integration points per dimension
    NC = NP^D             # total number of integration points
    X1, W1 = simplexquad(NP, collect(s)')
    @assert size(X1) == (NC, D)
    @assert size(W1) == (NC,)
    X = SMatrix{D, NC}(X1')
    dump(X)
    W = SVector{NC}(W1)

    setup = bernstein_setup(s)
    # XX = SVector{NC}([X[:,i] for i in 1:NC])
    # SVectors containing SVectors cause problems. We use the
    # "Identity" type to hide the inner SVectors while constructing
    # the outer ones.
    # XX = SVector{NC}(X[:,i] for i in 1:NC)
    XX = getindex.(SVector{NC}(Identity(X[:,i]) for i in 1:NC))
    ΛΛ = map(x -> cartesian2barycentric(setup, x), XX)
    Λ = SMatrix{N, NC}(ΛΛ[i][a] for a in 1:N, i in 1:NC)

    nα = binomial(P+D, D)
    bdb = Array{T}(undef, nα, nα)
    i2α = Array{SVector{N, Int}}(undef, nα)

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
            bdb[i,j] = integrate_λ(f, Λ, W)
        end
        @assert j == nα
    end
    @assert i == nα

    i2α, Symmetric(bdb)
end

bernstein_products(T, D, P) = bernstein_products(T, Val(D), Val(P))



@fastmath function integrate_λ(f::F, Λ::SMatrix{N, NC, T}, W::SVector{NC, T}
                               ) where {F, N, NC, T}
    @inbounds begin
        s = zero(W[1] * f(Λ[:, 1]))
        for n in 1:NC
            s += W[n] * f(Λ[:, n])
        end
        s
    end
end
