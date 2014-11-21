type BSplineBasis
    knots::Vector{Float64}
    order::UInt8

    function BSplineBasis(knots, order, extend=false)
        for (kn, kp) in zip(knots[2:end], knots[1:end-1])
            @assert kn >= kp "Knot vector must be nondecreasing"
        end
        if extend
            knots = [fill(knots[1], order-1), knots, fill(knots[end], order-1)]
        else
            for d in 1:order-1
                @assert(knots[1+d] == knots[1] && knots[end-d] == knots[end],
                        "Expected $order repeated knots on either end")
            end
        end
        return new(knots, order)
    end
end

dim(basis::BSplineBasis) = length(basis.knots) - basis.order


