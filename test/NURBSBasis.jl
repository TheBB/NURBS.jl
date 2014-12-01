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

# Single point, no derivatives
basis = NURBSBasis(BSplineBasis(0, 1, 2, 2))
testeval_snn(basis, 0, [1, 0], 1:2)
testeval_snn(basis, 0.5, [1, 0], 2:3)
testeval_snn(basis, 1, [0, 1], 2:3)

basis = NURBSBasis(BSplineBasis(0, 2, 2, 2), [1.5, 0.5, 1])
testeval_snn(basis, 0.2, [0.923076923077, 0.0769230769231], 1:2)
