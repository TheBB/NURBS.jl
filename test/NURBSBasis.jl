using NURBS.Utils
using NURBS.Bases


# NURBSBasis inner constructor
# ========================================================================

bs = BSplineBasis(linspace(0, 4, 5), 3)

basis = NURBSBasis(bs, [1.1, 1.2, 1.3, 1.4, 1.5, 1.6], 1)
@test basis.bs == bs
@test basis.weights == [1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
@test basis.deriv == 1

basis = NURBSBasis(bs)
@test basis.bs == bs
@test basis.weights == ones(6)
@test basis.deriv == 0

basis = NURBSBasis(bs, 2)
@test basis.bs == bs
@test basis.weights == ones(6)
@test basis.deriv == 2

basis = NURBSBasis(bs, linspace(4, 5, 6))
@test basis.bs == bs
@test basis.weights == linspace(4, 5, 6)
@test basis.deriv == 0

@test_throws ArgumentError NURBSBasis(deriv(bs))
@test_throws ArgumentError NURBSBasis(bs, linspace(1, 2, 5))
@test_throws ArgumentError NURBSBasis(bs, linspace(1, 2, 7))
@test_throws ArgumentError NURBSBasis(bs, [-1.0, ones(5)])
@test_throws ArgumentError NURBSBasis(bs, 3)
@test_throws ArgumentError NURBSBasis(BSplineBasis([0, 1], 2), 2)


# NURBSBasis length, size and domain
# ========================================================================

for bs in [BSplineBasis(0, 1, 4, 4), BSplineBasis(0, 2, 5, 3)]
    nbs = NURBSBasis(bs)
    @test length(nbs) == length(bs)
    @test size(nbs) == size(bs)
    @test domain(nbs) == domain(bs)
end


# NURBSBasis evaluation
# ========================================================================

t = sqrt(2)/2
wts = [1, t, 1, t, 1, t, 1, t, 1]
kts = [0, pi/2, pi/2, pi, pi, 3pi/2, 3pi/2, 2pi]
coeffs = [1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1; 1 0]
bs = NURBSBasis(BSplineBasis(kts, 3), wts)
norms = sum(bs(linspace(0, 2pi, 100), coeffs) .^ 2, 2)
@test_approx_eq norms ones(100)
