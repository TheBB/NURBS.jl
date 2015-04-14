import Iterators: product

using NURBS.Utils
using NURBS.Bases

info("Testing NURBS.Bases")


# BSplineBasis inner constructor
# ========================================================================

basis = BSplineBasis(linspace(0, 4, 5), 3)
@test basis.knots == [0, 0, 0, 1, 2, 3, 4, 4, 4]
@test basis.order == 3
@test basis.deriv == 0

basis = BSplineBasis([0, 0, 0, 1, 2, 3, 4, 4, 4], 3, 1; extend=false)
@test basis.knots == [0, 0, 0, 1, 2, 3, 4, 4, 4]
@test basis.order == 3
@test basis.deriv == 1

@test_throws ArgumentError BSplineBasis([2, 2, 1, 0, 0], 2, 0; extend=false)
@test_throws ArgumentError BSplineBasis([0, 0, 1, 2, 2], 3, 0; extend=false)
@test_throws ArgumentError BSplineBasis([0, 0, 1, 2, 3], 2, 0; extend=false)
@test_throws ArgumentError BSplineBasis([0, 1, 2, 3, 3], 2, 0; extend=false)
@test_throws ArgumentError BSplineBasis([0, 1], 2, 2)


# BSplineBasis outer constructors
# ========================================================================

basis = BSplineBasis(0, 1, 4, 4)
@test basis.knots == [0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1]
@test basis.order == 4
@test basis.deriv == 0


# BSplineBasis length, size and domain
# ========================================================================

basis = BSplineBasis(0, 1, 4, 4)
@test length(basis) == 7
@test size(basis) == (7,)
@test domain(basis) == Interval(0, 1)

basis = BSplineBasis(0, 2, 5, 3)
@test length(basis) == 7
@test size(basis) == (7,)
@test domain(basis) == Interval(0, 2)


# BSplineBasis evaluation
# ========================================================================

# Single point, no derivatives
function testeval_snn(b, x, vals, idxs)
    (rvals, ridxs) = b(x)
    @test ridxs == idxs
    @test_approx_eq rvals vals

    coeffs = rand(length(b))
    @test_approx_eq dot(coeffs[ridxs], vals) b(x, coeffs)

    for (v, i) in zip(vals, idxs)
        @test_approx_eq b[i](x) v
    end

    for i in 1:length(b)
        if i ∉ idxs
            @test_approx_eq b[i](x) 0.0
        end
    end
end

basis = BSplineBasis(0, 1, 2, 2)
testeval_snn(basis, 0, [1, 0], 1:2)
testeval_snn(basis, 0.2, [0.6, 0.4], 1:2)
testeval_snn(basis, 0.4, [0.2, 0.8], 1:2)
testeval_snn(basis, 0.5, [1, 0], 2:3)
testeval_snn(basis, 0.75, [0.5, 0.5], 2:3)
testeval_snn(basis, 1, [0, 1], 2:3)

basis = BSplineBasis(0, 1, 1, 3)
testeval_snn(basis, 0, [1, 0, 0], 1:3)
testeval_snn(basis, 0.3, [0.49, 0.42, 0.09], 1:3)
testeval_snn(basis, 0.8, [0.04, 0.32, 0.64], 1:3)

basis = BSplineBasis(0, 2, 2, 4)
testeval_snn(basis, 0.2, [0.512, 0.434, 0.052, 0.002], 1:4)
testeval_snn(basis, 0.7, [0.027, 0.49525, 0.392, 0.08575], 1:4)
testeval_snn(basis, 1.0, [0.25, 0.5, 0.25, 0], 2:5)
testeval_snn(basis, 1.5, [0.03125, 0.25, 0.59375, 0.125], 2:5)
testeval_snn(basis, 2.0, [0, 0, 0, 1], 2:5)

basis = BSplineBasis([0, 1, 3], 4)
testeval_snn(basis, 0.5, [0.125, 0.680555555556, 0.18055555555555564, 0.013888888888888902], 1:4)
testeval_snn(basis, 1.4, [0.227555555555557, 0.483555555556, 0.280888888889, 0.008], 2:5)
testeval_snn(basis, 2.1, [0.0405, 0.26325, 0.529875, 0.166375], 2:5)


# Multiple points, no derivatives
function testeval_mnn(b, xs)
    res = b(xs)
    for (x, (vals, idxs)) in zip(xs, res)
        @test typeof(vals) <: Vector
        testeval_snn(b, x, vals, idxs)
    end

    coeffs = rand(length(b), 3)
    resc = b(xs, coeffs)
    for (i, (vals, idxs)) in enumerate(res)
        @test_approx_eq vals'*coeffs[idxs,:] resc[i,:]
    end
end

basis = BSplineBasis(0, 2, 2, 2)
pts = [0.0, 0.4, 0.8, 1.0, 1.5, 2.0]
testeval_mnn(basis, pts)
@test_approx_eq basis[1](pts) [1.0; 0.6; 0.2; zeros(3)]
@test_approx_eq basis[2](pts) [0.0, 0.4, 0.8, 1.0, 0.5, 0.0]
@test_approx_eq basis[3](pts) [zeros(4); 0.5; 1.0]

basis = BSplineBasis(0, 3, 4, 3)
pts = [0.0, 0.4, 0.8, 1.0, 1.5, 2.0, 2.2, 2.8, 3.0]
testeval_mnn(basis, pts)
@test_approx_eq basis[1](pts) [1.0; 0.2177777777777774; zeros(7)]
@test_approx_eq basis[3](pts) [0.0, 0.1422222222222222, 0.56222222222222222, 0.722222222222222,
                               0.5, 0.0555555555555555, 0.00222222222222222, 0.0, 0.0]
@test_approx_eq basis[6](pts) [zeros(7); 0.53777777777777777; 1.0]

basis = BSplineBasis([0, 1, 3, 4], 4)
pts = [0.0, 0.8, 1.6, 2.4, 3.2, 4.0]
testeval_mnn(basis, pts)
@test_approx_eq basis[2](pts) [0.0; 0.5795555555555555; 0.15244444444444444; 0.012; zeros(2)]
@test_approx_eq basis[4](pts) [0.0, 0.0426666666666666, 0.29333333333333333, 0.542222222222222,
                               0.3697777777777777, 0.0]


# Single point, derivatives
function testeval_sdn(b, x, dt=1e-6, tol=1e-8)
    dvals, didxs = deriv(b)(x)
    (dvp, didp), (dvn, didn) = b([x+dt, x-dt])
    @test didxs == didp == didn
    nderiv = (dvp - dvn) / 2dt
    @test_approx_eq_eps nderiv dvals tol

    coeffs = rand(length(b), 4)
    @test_approx_eq dvals'*coeffs[didxs,:] deriv(b)(x, coeffs)

    for (v, i) in zip(dvals, didxs)
        @test_approx_eq deriv(b[i])(x) v
    end

    for i in 1:length(b)
        if !(i in didxs)
            @test_approx_eq deriv(b[i])(x) 0.0
        end
    end
end

basis = BSplineBasis(0, 1, 2, 2)
for pt in [0.25, 0.51, 0.9]
    testeval_sdn(basis, pt)
end
basis = deriv(basis)
@test_throws ArgumentError deriv(basis)

basis = BSplineBasis([0, 0.5, 0.7, 0.8], 3)
for i in 1:2
    for pt in [0.1, 0.3, 0.6, 0.73]
        testeval_sdn(basis, pt)
    end
    basis = deriv(basis)
end
@test_throws ArgumentError deriv(basis)

basis = BSplineBasis(0, 2, 4, 4)
for i in 1:3
    for pt in [0.1, 0.7, 1.2, 1.6]
        testeval_sdn(basis, pt)
    end
    basis = deriv(basis)
end
@test_throws ArgumentError deriv(basis)


# Multiple points, derivatives
function testeval_mdn(b, xs, dt=1e-6, tol=1e-8)
    res = deriv(b)(xs)
    resn = b(xs+dt)
    resp = b(xs-dt)
    for ((v, id), (nv, nid), (pv, pid)) in zip(res, resn, resp)
        @test typeof(v) <: Vector
        @test id == nid == pid
        nderiv = (nv - pv) / 2dt
        @test_approx_eq_eps nderiv v tol
    end

    coeffs = rand(length(b))
    resc = deriv(b)(xs, coeffs)
    for (i, (v, id)) in enumerate(res)
        @test_approx_eq dot(v, coeffs[id]) resc[i]
    end
end

basis = BSplineBasis([0, 0.1, 0.8, 1.2], 2)
testeval_mdn(basis, [0.02, 0.3, 0.7, 0.9])

basis = BSplineBasis(6, 12, 3, 3)
for i in 1:2
    testeval_mdn(basis, [6.43, 7.8, 8.3, 9.1, 10.1, 11.8])
    basis = deriv(basis)
end

basis = BSplineBasis(5, 15, 10, 4)
for i in 1:3
    testeval_mdn(basis, collect(linspace(5.5, 14.5, 10)))
    basis = deriv(basis)
end


# BSplineBasis indexing
# ========================================================================

basis = BSplineBasis(0, 4, 4, 2)
@test_throws BoundsError basis[0]
@test_throws BoundsError basis[6]
@test basis[1].basis == basis
@test basis[2].index == 2
@test basis[5].deriv == 0
@test domain(basis[1]) == Interval(0, 1)
@test domain(basis[2]) == Interval(0, 2)
@test domain(basis[5]) == Interval(3, 4)
@test deriv(basis[1]).basis == basis
@test deriv(basis[2]).index == 2
@test deriv(basis[3]).deriv == 1
@test_throws ArgumentError deriv(deriv(basis[2]))

basis = BSplineBasis([0, 4, 5, 6], 4)
@test_throws BoundsError basis[7]
@test basis[6].basis == basis
@test basis[2].deriv == 0
@test basis[1].index == 1
@test domain(basis[1]) == Interval(0, 4)
@test domain(basis[2]) == Interval(0, 5)
@test domain(basis[3]) == Interval(0, 6)
@test domain(basis[4]) == Interval(0, 6)
@test domain(basis[5]) == Interval(4, 6)
@test domain(basis[6]) == Interval(5, 6)
@test deriv(basis[1]).basis == basis
@test deriv(basis[4]).index == 4
@test deriv(deriv(basis[4])).deriv == 2
@test_throws ArgumentError deriv(deriv(deriv(deriv(basis[3]))))

basis = deriv(basis)
@test deriv(basis[2]).deriv == 2
@test deriv(basis[2]).basis.deriv == 1
@test_throws ArgumentError deriv(deriv(deriv(basis[1])))
