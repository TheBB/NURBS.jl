immutable NURBSBasis <: Basis1D
    bs::BSplineBasis
    weights::Vector{Float64}
    deriv::Int

    function NURBSBasis(bs::BSplineBasis, weights, deriv=0)
        @assert_ex(bs.deriv == 0,
                   ArgumentError("Underlying B-Spline basis must not be differentiated"))
        @assert_ex(length(bs) == length(weights),
                   ArgumentError("Number of weights must be equal to number of basis functions"))
        @assert_ex(minimum(weights) > 0.0, ArgumentError("Weights must be positive"))
        @assert_ex(deriv <= min(nderivs(bs), 2), ArgumentError("Differentiation order not supported"))
        new(bs, weights, deriv)
    end

    NURBSBasis(b::BSplineBasis, deriv::Int) = NURBSBasis(b, ones(length(b)), deriv)
    NURBSBasis(b::BSplineBasis) = NURBSBasis(b, ones(length(b)))
end

typealias NURBS BasisFunction1D{NURBSBasis}


# Inherit functions from B-Spline bases
Base.length(b::NURBSBasis) = length(b.bs)
Base.size(b::NURBSBasis) = size(b.bs)
Base.getindex(b::NURBSBasis, i) = NURBS(b, i, b.deriv)

nderivs(b::NURBSBasis) = min(2, nderivs(b.bs))

domain(b::NURBSBasis) = domain(b.bs)
domain(b::NURBS) = domain(b.basis[b.index])

order(b::NURBSBasis) = throw(ArgumentError("NURBS are not polynomial"))
order(b::NURBS) = throw(ArgumentError("NURBS are not polynomial"))

deriv(b::NURBSBasis) = NURBSBasis(b.bs, b.weights, b.deriv + 1)
deriv(b::NURBSBasis, order) = NURBSBasis(b.bs, b.weights, b.deriv + order)

supported(b::NURBSBasis, pts) = supported(b.bs, pts)

function evaluate_raw{T<:Real}(b::NURBSBasis, pts::Vector{T}, deriv::Int, rng::UnitRange{Int})
    const bwts = b.weights[rng]

    bvals = evaluate_raw(b.bs, pts, 0, rng) .* bwts
    wts = sum(bvals, 1)

    if deriv == 0
        return bvals ./ wts
    end

    bvals1 = evaluate_raw(b.bs, pts, 1, rng) .* bwts
    wts1 = sum(bvals1, 1)
    d1 = (bvals1 .* wts - bvals .* wts1) ./ (wts .^ 2)

    if deriv == 1
        return d1
    end

    bvals2 = evaluate_raw(b.bs, pts, 2, rng) .* bwts
    wts2 = sum(bvals2, 1)
    d2 = (bvals2 .* wts - bvals .* wts2) ./ (wts .^ 2)

    if deriv == 2
        return d2 - 2 * d1 .* wts1 ./ wts
    end

    throw(ArgumentError("Third order derivatives or higher are not supported by NURBS bases"))
end
