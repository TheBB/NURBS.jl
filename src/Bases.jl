module Bases

import Iterators: groupby, imap

export Basis, dim, supported, evaluate_raw, evaluate,
       BSplineBasis, uniform_bsbasis


# Abstract Basis type
# ========================================================================
#
# A subtype MyBasis <: Basis should implement the following functions.
#
# dim(b::MyBasis)
#     Dimension of the basis.
#
# supported(b::MyBasis, pts::Vector{Float64})
#     Organizes the points into subsets on which the same basis functions are supported.
#     Return any iterator of (Vector{Float64}, Range).
#     May assume that pts is sorted.
#
# evaluate_raw(b::MyBasis, pt::Float64, rng::Range)
#     Evaluate the basis functions in rng at pt.
#     May assume that rng == supported(b, pt).

abstract Basis

evaluate(b::Basis, pt::Float64) = evaluate_raw(b, pt, supported(b, pt))
supported(b::Basis, pt::Float64) = first(supported(b, [pt]))[2]


# B-Spline basis
# ========================================================================

type BSplineBasis <: Basis
    knots::Vector{Float64}
    order::Int

    function BSplineBasis(knots, order, extend=true)
        for (kn, kp) in zip(knots[2:end], knots[1:end-1])
            @assert(kn >= kp, "Knot vector must be nondecreasing")
        end

        if extend
            knots = [fill(knots[1], order-1), knots, fill(knots[end], order-1)]
        else
            for d in 1:order-1
                @assert(knots[1+d] == knots[1] && knots[end-d] == knots[end],
                        "Expected $order repeated knots on either end")
            end
        end

        new(knots, order)
    end
end

uniform_bsbasis(elements, order=4, lower=0.0, upper=1.0) =
    BSplineBasis(linspace(lower, upper, elements+1), order)

dim(basis::BSplineBasis) = length(basis.knots) - basis.order

function supported(b::BSplineBasis, pts::Vector{Float64})
    @assert(b.knots[1] <= pts[1] <= pts[end] <= b.knots[end],
            "Evaluation points must lie within parameter domain")

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

function evaluate_raw(b::BSplineBasis, pt::Float64, rng::UnitRange{Int})
    # Basis values of order 1 (piecewise constants)
    bvals = zeros(Float64, b.order)
    bvals[end] = 1.0

    const p = b.order
    const bi = rng.start + p

    # Increase order
    for k in 0:p-2
        dxs = b.knots[bi : bi+k] - b.knots[bi-k-1 : bi-1]
        bvals[p-k:end] ./= dxs

        lft = pt - b.knots[bi-k-2:bi-1]
        rgt = b.knots[bi:bi+k] - pt
        bvals[p-k-1:end] =
            lft .* bvals[p-k-1:end] + [rgt .* bvals[p-k:end], 0]
    end

    bvals
end

end  # module BSplineBasis
