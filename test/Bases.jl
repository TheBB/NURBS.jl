import Iterators: product

using NURBS.Bases

info("Testing NURBS.Bases")


# BSplineBasis inner constructor
# ========================================================================

basis = BSplineBasis(linspace(0, 4, 5), 3)
@test(basis.order == 3)
@test(basis.knots == [0, 0, 0, 1, 2, 3, 4, 4, 4])
@test(dim(basis) == 6)

basis = BSplineBasis([0, 0, 0, 1, 2, 3, 4, 4, 4], 3, false)
@test(basis.order == 3)
@test(basis.knots == [0, 0, 0, 1, 2, 3, 4, 4, 4])
@test(dim(basis) == 6)

@test_throws(ErrorException, BSplineBasis([2, 2, 1, 0, 0], 2, false))
@test_throws(ErrorException, BSplineBasis([0, 0, 1, 2, 2], 3, false))
@test_throws(ErrorException, BSplineBasis([0, 0, 1, 2, 3], 2, false))
@test_throws(ErrorException, BSplineBasis([0, 1, 2, 3, 3], 2, false))


# BSplineBasis outer constructors
# ========================================================================

basis = uniform_bsbasis(4, 4, 0, 1)
@test(basis.order == 4)
@test(basis.knots == [0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1])
@test(dim(basis) == 7)


# BSplineBasis evaluation
# ========================================================================

basis = uniform_bsbasis(2, 2, 0, 1)
@test_approx_eq(evaluate(basis, 0.0), [1, 0])
@test_approx_eq(evaluate(basis, 0.1), [0.8, 0.2])
@test_approx_eq(evaluate(basis, 0.4), [0.2, 0.8])
@test_approx_eq(evaluate(basis, 0.5), [1, 0])
@test_approx_eq(evaluate(basis, 0.75), [0.5, 0.5])

basis = uniform_bsbasis(2, 2, 0, 1)
@test(supported(basis, 0.0) == 1:2)
@test(supported(basis, 0.2) == 1:2)
@test(supported(basis, 0.5) == 2:3)
@test(supported(basis, 0.8) == 2:3)
@test(supported(basis, 1.0) == 2:3)

basis = uniform_bsbasis(1, 3, 0, 1)
@test_approx_eq(evaluate(basis, 0.0), [1, 0, 0])
@test_approx_eq(evaluate(basis, 0.3), [0.49, 0.42, 0.09])
@test_approx_eq(evaluate(basis, 0.8), [0.04, 0.32, 0.64])

basis = uniform_bsbasis(2, 3, 0, 1)
@test(supported(basis, 0.0) == 1:3)
@test(supported(basis, 0.4) == 1:3)
@test(supported(basis, 0.5) == 2:4)
@test(supported(basis, 1.0) == 2:4)

basis = uniform_bsbasis(1, 4, 0, 2)
@test_approx_eq(evaluate(basis, 0.0), [1, 0, 0, 0])
@test_approx_eq(evaluate(basis, 0.7), [0.274625, 0.443625, 0.238875, 0.042875])
@test_approx_eq(evaluate(basis, 1.4), [0.027, 0.189, 0.441, 0.343])

basis = uniform_bsbasis(2, 4, 0, 2)
@test(supported(basis, 0.0) == 1:4)
@test(supported(basis, 1.0) == 2:5)
@test(supported(basis, 1.9) == 2:5)
@test(supported(basis, 2.0) == 2:5)

function test_derivatives(basis, pt, dt=1e-6, tol=1e-8)
    lft = evaluate(basis, pt - dt, basis.order-1)
    rgt = evaluate(basis, pt + dt, basis.order-1)
    mid = evaluate(basis, pt, basis.order-1)
    vals = (rgt - lft) / 2dt
    @test_approx_eq_eps(mid[:,2:end], vals[:,1:end-1], tol)
end

for (ord, pt) in product(2:4, [0.25, 0.47, 0.89])
    test_derivatives(uniform_bsbasis(2, ord, 0, 1), pt)
end
