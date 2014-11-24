type BSplineBasis{T}
    knots::Vector{T}
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

dim(knots::BSplineBasis) = length(basis.knots) - basis.order

function evaluate_raw{T}(knots::Vector{T}, coeffs::Vector{T}, point::T)
    n = length(coeffs)
    while n > 1
        α = (point - knots[1:n-1]) ./ (knots[end-n+2:end] - knots[1:n-1])
        coeffs = (1 - α) .* coeffs[1:end-1] + α .* coeffs[2:end]
        knots = knots[2:end-1]
        n -= 1
    end

    return coeffs
end
