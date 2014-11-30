module Bases

import Iterators: groupby, imap
import ..Utils: Interval

export Basis, Basis1D, domain, deriv, order,
       BSplineBasis


# Abstract Basis types
# ========================================================================

abstract Basis
abstract Basis1D <: Basis


# B-Spline basis
# ========================================================================

type BSplineBasis <: Basis1D
    knots::Vector{Float64}
    order::Int
    deriv::Int

    function BSplineBasis(knots, order, deriv=0, extend=true)
        @assert(deriv < order)

        for (kn, kp) in zip(knots[2:end], knots[1:end-1])
            @assert(kn >= kp)
        end

        if extend
            knots = [fill(knots[1], order-1), knots, fill(knots[end], order-1)]
        else
            for d in 1:order-1
                @assert(knots[1+d] == knots[1] && knots[end-d] == knots[end])
            end
        end

        new(knots, order, deriv)
    end
end

BSplineBasis(lft::Real, rgt::Real, elements::Int, order::Int) =
    BSplineBasis(linspace(lft, rgt, elements+1), order, 0)

type BSpline
    basis::BSplineBasis
    index::Int
    deriv::Int

    function BSpline(basis, index, deriv)
        @assert(1 <= index <= length(basis))
        @assert(0 <= deriv < basis.order - basis.deriv)
        new(basis, index, deriv)
    end

    BSpline(basis, index) = BSpline(basis, index, 0)
end

Base.length(b::BSplineBasis) = length(b.knots) - b.order
Base.size(b::BSplineBasis) = (length(b),)
Base.getindex(b::BSplineBasis, i) = BSpline(b, i, 0)

domain(b::BSplineBasis) = Interval(b.knots[1], b.knots[end])
domain(b::BSpline) = Interval(b.basis.knots[b.index], b.basis.knots[b.index+b.basis.order])

order(b::BSplineBasis) = b.order - b.deriv
order(b::BSpline) = b.basis.order - b.basis.deriv - b.deriv

deriv(b::BSplineBasis) = BSplineBasis(b.knots, b.order, b.deriv+1, false)
deriv(b::BSplineBasis, order) = BSplineBasis(b.knots, b.order, b.deriv+order, false)
deriv(b::BSpline) = BSpline(b.basis, b.index, b.deriv+1)
deriv(b::BSpline, order) = BSpline(b.basis, b.index, b.deriv+order)

function Base.call{T<:Real}(b::BSplineBasis, pt::T)
    rng = supported(b, pt)
    (squeeze(evaluate_raw(b, [pt], b.deriv, rng), 2), rng)
end

function Base.call{T<:Real}(b::BSpline, pt::T)
    rng = supported(b.basis, pt)
    if b.index ∉ rng return 0.0 end
    evaluate_raw(b.basis, [pt], b.deriv + b.basis.deriv, rng)[1 + b.index - rng.start, 1]
end

function Base.call{T<:Real}(b::BSplineBasis, pts::Vector{T})
    res = (Vector{Float64}, UnitRange{Int})[]
    sizehint(res, length(pts))

    for (subpts, rng) in supported(b, pts)
        out = evaluate_raw(b, subpts, b.deriv, rng)
        for i in 1:length(subpts)
            push!(res, (out[:,i], rng))
        end
    end

    res
end

function Base.call{T<:Real}(b::BSpline, pts::Vector{T})
    res = zeros(Float64, length(pts))

    i = 1
    for (subpts, rng) in supported(b.basis, pts)
        i += length(subpts)
        if b.index ∉ rng continue end

        out = evaluate_raw(b.basis, subpts, b.basis.deriv + b.deriv, rng)
        res[i-length(subpts):i-1] = out[findin(rng, b.index), :]
    end

    res
end

function Base.call{S<:Real, T<:Real}(b::BSplineBasis, pt::S, coeffs::Vector{T})
    (vals, idxs) = b(pt)
    dot(vals, coeffs[idxs])
end

function Base.call{S<:Real, T<:Real}(b::BSplineBasis, pt::S, coeffs::Matrix{T})
    (vals, idxs) = b(pt)
    vals' * coeffs[idxs,:]
end

Base.call{S<:Real, T<:Real}(b::BSplineBasis, pts::Vector{S}, coeffs::Vector{T}) =
    Float64[dot(vals, coeffs[idxs]) for (vals, idxs) in b(pts)]

function Base.call{S<:Real, T<:Real}(b::BSplineBasis, pts::Vector{S}, coeffs::Matrix{T})
    res = zeros(Float64, length(pts), size(coeffs, 2))

    for (i, (vals, idxs)) in enumerate(b(pts))
        res[i,:] = vals' * coeffs[idxs,:]
    end

    res
end

function supported{T<:Real}(b::BSplineBasis, pt::T)
    kidx = b.order - 1 + searchsorted(b.knots[b.order:end], pt).stop
    stop = b.knots[kidx] == b.knots[end] ? kidx - b.order : kidx
    stop - b.order + 1 : stop
end

function supported{T<:Real}(b::BSplineBasis, pts::Vector{T})
    (min, max) = extrema(pts)
    @assert(min in domain(b) && max in domain(b))

    idxs = zeros(Int, length(pts))

    if !issorted(pts)
        for (i, pt) in enumerate(pts)
            idxs[i] = supported(b, pt).stop
        end
    else
        kidx = b.order
        for (i, pt) in enumerate(pts)
            kidx = kidx - 1 + searchsorted(b.knots[kidx:end], pt).stop
            idxs[i] = b.knots[kidx] == b.knots[end] ? kidx - b.order : kidx
        end
    end

    imap(groupby(enumerate(idxs), i -> i[2])) do i
        (pts[i[1][1]:i[end][1]], i[1][2] - b.order + 1 : i[1][2])
    end
end

macro bs_er_scale(bvals, knots, mid, num)
    :($bvals ./= $knots[$mid:$mid+$num] - $knots[$mid-$num-1:$mid-1])
end

function evaluate_raw{T<:Real}(b::BSplineBasis, pts::Vector{T}, deriv::Int, rng::UnitRange{Int})
    # Basis values of order 1 (piecewise constants)
    bvals = zeros(Float64, (b.order, length(pts)))
    bvals[end,:] = 1.0

    const p = b.order
    const bi = rng.start + p

    # Order increment
    for k in 0:p-deriv-2
        @bs_er_scale(bvals[p-k:end,:], b.knots, bi, k)

        for (i, kp, kn) in zip(p-k-1:p-1, b.knots[bi-k-2:bi-2], b.knots[bi:bi+k])
            bvals[i,:] .*= (pts - kp)'
            bvals[i,:] += bvals[i+1,:] .* (kn - pts)'
        end
        bvals[end,:] .*= (pts - b.knots[bi-1])'
    end

    # Differentiation
    for k = p-deriv-1:p-2
        @bs_er_scale(bvals[p-k:end,:], b.knots, bi, k)

        bvals[1:end-1,:] = -diff(bvals, 1)
        bvals *= k + 1
    end

    bvals
end

end  # module Bases
