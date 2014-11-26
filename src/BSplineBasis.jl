module Basis

export BSplineBasis, uniform_basis, dim, evaluate_raw

type BSplineBasis
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

uniform_basis(elements, order=4, lower=0.0, upper=1.0) =
    BSplineBasis(linspace(lower, upper, elements+1), order)

dim(basis::BSplineBasis) = length(basis.knots) - basis.order

function evaluate_raw(basis::BSplineBasis, point::Float64)
    knots, order = basis.knots, basis.order
    
    # Point should satisfy knots[stop] <= point < knots[stop+1],
    # unless at the right hand side of the interval
    stop = searchsorted(knots, point).stop
    if knots[stop] == knots[end]
        stop -= order
    end
    start = stop - order + 1

    # Basis values of order 1 (piecewise constants)
    bvals = zeros(Float64, order)
    bvals[end] = 1.0

    for k in 2:order
        for i in order-k+2:order
            dx = knots[start+i+k-2] - knots[start+i-1]
            if dx > 0
                bvals[i] /= dx
            else
                bvals[i] = 0.0
            end
        end

        for i in order-k+1:order
            bvals[i] *= (point - knots[start+i-1])
            if i < order
                bvals[i] += (knots[start+i+k-1] - point) * bvals[i+1]
            end
        end
    end

    return bvals, start:stop
end

end  # module BSplineBasis
