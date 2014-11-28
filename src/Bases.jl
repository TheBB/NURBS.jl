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
        @assert(0 <= deriv < basis.order)
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

function Base.call(b::BSplineBasis, pt::Real)
    rng = supported(b, pt)
    (squeeze(evaluate_raw(b, pt, b.deriv, rng)[:, 1 + b.deriv], 2), rng)
end

function Base.call(b::BSpline, pt::Real)
    rng = supported(b.basis, pt)
    if b.index ∉ rng return 0.0 end
    der = b.deriv + b.basis.deriv
    evaluate_raw(b.basis, pt, der, rng)[1 + b.index - rng.start, 1 + der]
end

Base.call(b::BSplineBasis, pts) = [b(pt) for pt in pts]
Base.call(b::BSpline, pts) = [b(pt) for pt in pts]

Base.call(b::BSplineBasis, pts, coeffs::Vector) =
    [dot(coeffs[idxs], vals) for (vals, idxs) in b(pts)]

Base.call(b::BSplineBasis, pts, coeffs) =
    [vals' * coeffs[idxs,:] for (vals, idxs) in b(pts)]

function supported(b::BSplineBasis, pts::Vector{Float64})
    @assert(b.knots[1] <= pts[1] <= pts[end] <= b.knots[end])

    idxs = zeros(Int, length(pts))

    kidx = b.order
    for (i, pt) in enumerate(pts)
        kidx = kidx - 1 + searchsorted(b.knots[kidx:end], pt).stop
        idxs[i] = b.knots[kidx] == b.knots[end] ? kidx - b.order : kidx
    end

    imap(groupby(enumerate(idxs), i -> i[2])) do i
        (pts[i[1][1]:i[end][1]], i[1][2] - b.order + 1 : i[1][2])
    end
end

supported(b::BSplineBasis, pt::Real) = first(supported(b, Float64[pt]))[2]

function evaluate_raw(b::BSplineBasis, pt::Real, nder::Int, rng::UnitRange{Int})
    # Basis values of order 1 (piecewise constants)
    bvals = zeros(Float64, (b.order,  nder+1))
    bvals[end,end] = 1.0

    const p = b.order
    const bi = rng.start + p
    col = nder + 1

    # Iterate over orders
    for k in 0:p-2
        # Scale basis functions
        dxs = b.knots[bi : bi+k] - b.knots[bi-k-1 : bi-1]
        bvals[p-k:end, col:end] ./= dxs

        if k > p-2-nder
            # Copy scaled basis functions to next level
            bvals[:, col-1] = bvals[:, col]

            # Differentiate the remainder
            # This 'unscales' the functions
            bvals[:, col:end] = [bvals[1:end-1, col:end] - bvals[2:end, col:end];
                                 bvals[end, col:end]] * (k + 1)

            col -= 1
        end

        # Apply order increment formula to the highest level
        # This also 'unscales' the functions
        lft = (pt - b.knots[bi-k-2:bi-1]) .* bvals[p-k-1:end, col]
        rgt = [(b.knots[bi:bi+k] - pt) .* bvals[p-k:end, col], 0]
        bvals[p-k-1:end, col] = lft + rgt
    end

    bvals
end

end  # module Bases
