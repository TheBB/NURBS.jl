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
#     Returns any iterator of (Vector{Float64}, Range).
#     May assume that pts is sorted.
#
# evaluate_raw(b::MyBasis, pt::Float64, nder::Int, rng::Range)
#     Evaluate the basis functions in rng at pt, as well as their first nder derivatives.
#     Returns an array of dimensions |rng| × (nder+1).
#     May assume that rng == supported(b, pt).

abstract Basis

evaluate(b::Basis, pt::Float64, nder::Int) = evaluate_raw(b, pt, nder, supported(b, pt))
evaluate(b::Basis, pt::Float64) = evaluate_raw(b, pt, 0, supported(b, pt))[:]
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

function evaluate_raw(b::BSplineBasis, pt::Float64, nder::Int, rng::UnitRange{Int})
    @assert(nder < b.order, "Higher order derivatives not supported")

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
