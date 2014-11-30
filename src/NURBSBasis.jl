immutable NURBSBasis <: Basis1D
    bs::BSplineBasis
    weights::Vector{Float64}
    deriv::Int

    function NURBSBasis(basis::BSplineBasis, weights, deriv)
        @assert(bs.deriv == 0)
        @assert(length(basis) == length(weights))
        @assert(deriv <= nderivs(basis))
        new(basis, weights, deriv)
    end

    NURBSBasis{T<:Real}(b::BSplineBasis, weights::Vector{T}) = NURBSBasis(b, weights, 0)
    NURBSBasis(b::BSplineBasis, deriv::Int) = NURBSBasis(b, ones(length(b)), deriv)
    NURBSBasis(b::BSplineBasis) = NURBSBasis(b, ones(length(b)), 0)
end

typealias NURBS BasisFunction1D{NURBSBasis}


# Inherit functions from B-Spline bases
Base.length(b::NURBSBasis) = length(b.bs)
Base.size(b::NURBSBasis) = size(b.bs)
Base.getindex(b::NURBSBasis) = NURBS(b, i, 0)

nderivs(b::NURBSBasis) = nderivs(b.bs) - b.deriv

domain(b::NURBSBasis) = domain(b.bs)
domain(b::NURBS) = domain(b.basis[b.index])

order(b::NURBSBasis) = throw(ArgumentError("NURBS are not polynomial"))
order(b::NURBS) = throw(ArgumentError("NURBS are not polynomial"))

deriv(b::NURBSBasis) = NURBSBasis(b.bs, b.weights, b.deriv + 1)
deriv(b::NURBSBasis, order) = NURBSBasis(b.bs, b.weights, b.deriv + order)

supported{T<:Real}(b::NURBSBasis, pt::T) = supported(b.bs, pt)
supported{T<:Real}(b::NURBSBasis, pts::Vector{T}) = supported(b.bs, pts)

function evaluate_raw{T<:Real}(b::NURBSBasis, pts::Vector{T}, deriv::Int, rng::UnitRange{Int})
    bvals = evaluate_raw(b.bs, pts, 0, rng)
    wts = b.weights[rng]' * bvals

    if deriv == 0
        return bvals .* b.weights[rng] ./ wts
    end

    bvals1 = evaluate_raw(b.bs, pts, 1, rng)
    wts1 = b.weights[rng]' * bvals
    d1 = (bvals1 .* wts - bvals .* wts1) ./ (wts .^ 2)

    if deriv == 1
        return d1
    end

    bvals2 = evaluate_raw(b.bs, pts, 2, rng)
    wts2 = b.weights[rng]' * bvals
    d2 = (bvals2 .* wts - bvals .* wts2) ./ (wts .^ 2)

    if deriv == 2
        return d2 - 2 * d1 .* wts1 ./ wts
    end

    throw(ArgumentError("Third order derivatives or higher are not supported by NURBS bases"))
end
