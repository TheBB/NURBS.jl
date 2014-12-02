module Bases

import Iterators: groupby, imap
import ..Utils: Interval

export Basis, Basis1D, domain, deriv, order,
       BSplineBasis,
       NURBSBasis


macro assert_ex(condition, error)
    :($condition ? nothing : throw($error))
end


# Abstract Basis types
# ========================================================================

abstract Basis
abstract Basis1D <: Basis

type BasisFunction1D{B<:Basis1D}
    basis::B
    index::Int
    deriv::Int

    function BasisFunction1D(basis, index, deriv)
        @assert_ex(1 <= index <= length(basis), BoundsError())
        @assert_ex(0 <= deriv <= nderivs(basis),
                   ArgumentError("Differentiation order not supported"))
        new(basis, index, deriv)
    end

    BasisFunction1D(basis, index) = BasisFunction1D(basis, index, 0)
end

deriv(b::BasisFunction1D) = typeof(b)(b.basis, b.index, b.deriv+1)
deriv(b::BasisFunction1D, order) = typeof(b)(b.basis, b.index, b.deriv+order)

function Base.call{B<:Basis1D, T<:Real}(b::B, pt::T)
    rng = supported(b, pt)
    (squeeze(evaluate_raw(b, [pt], b.deriv, rng), 2), rng)
end

function Base.call{B<:BasisFunction1D, T<:Real}(b::B, pt::T)
    rng = supported(b.basis, pt)
    if b.index ∉ rng return 0.0 end
    evaluate_raw(b.basis, [pt], b.deriv, rng)[1 + b.index - rng.start, 1]
end

function Base.call{B<:Basis1D, T<:Real}(b::B, pts::Vector{T})
    tp = (Vector{Float64}, UnitRange{Int})
    res = Array(tp, length(pts))

    j = 0
    for (subpts, rng) in supported(b, pts)
        out = evaluate_raw(b, subpts, b.deriv, rng)
        for i in 1:length(subpts)
            res[j+=1] = (out[:,i], rng)
        end
    end

    res
end

function Base.call{B<:BasisFunction1D, T<:Real}(b::B, pts::Vector{T})
    res = zeros(Float64, length(pts))

    i = 1
    for (subpts, rng) in supported(b.basis, pts)
        i += length(subpts)
        if b.index ∉ rng continue end

        out = evaluate_raw(b.basis, subpts, b.deriv, rng)
        res[i-length(subpts):i-1] = out[findin(rng, b.index), :]
    end

    res
end

function Base.call{B<:Basis1D, S<:Real, T<:Real}(b::B, pt::S, coeffs::Vector{T})
    (vals, idxs) = b(pt)
    dot(vals, coeffs[idxs])
end

function Base.call{B<:Basis1D, S<:Real, T<:Real}(b::B, pt::S, coeffs::Matrix{T})
    (vals, idxs) = b(pt)
    vals' * coeffs[idxs,:]
end

Base.call{B<:Basis1D, S<:Real, T<:Real}(b::B, pts::Vector{S}, coeffs::Vector{T}) =
    Float64[dot(vals, coeffs[idxs]) for (vals, idxs) in b(pts)]

function Base.call{B<:Basis1D, S<:Real, T<:Real}(b::B, pts::Vector{S}, coeffs::Matrix{T})
    res = zeros(Float64, length(pts), size(coeffs, 2))

    for (i, (vals, idxs)) in enumerate(b(pts))
        res[i,:] = vals' * coeffs[idxs,:]
    end

    res
end


# Include specific bases
# ========================================================================

include("BSplineBasis.jl")
include("NURBSBasis.jl")


end  # module Bases
