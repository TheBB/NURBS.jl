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
@test basis.weights == collect(linspace(4, 5, 6))
@test basis.deriv == 0

@test_throws ArgumentError NURBSBasis(deriv(bs))
@test_throws ArgumentError NURBSBasis(bs, linspace(1, 2, 5))
@test_throws ArgumentError NURBSBasis(bs, linspace(1, 2, 7))
@test_throws ArgumentError NURBSBasis(bs, [-1.0; ones(5)])
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

testeval_snn(bs, 0, [1, 0, 0], 1:3)
testeval_snn(bs, 0.5, [0.53237051599, 0.351554065456, 0.116075418554], 1:3)
testeval_snn(bs, 1.0, [0.1527439727714, 0.378440763761, 0.468815263468], 1:3)
testeval_snn(bs, 1.5, [0.0020838736033, 0.062440576626, 0.935475549771], 1:3)
testeval_snn(bs, 2.0, [0.597709490779, 0.317802576968, 0.0844879322529], 3:5)
testeval_snn(bs, 2.5, [0.1943378565224, 0.398037455731, 0.407624687747], 3:5)
testeval_snn(bs, 3.0, [0.00853541363591, 0.1218408648272, 0.869623721537], 3:5)
testeval_snn(bs, 3.5, [0.664247692256, 0.277702618453, 0.0580496892913], 5:7)
testeval_snn(bs, 4.0, [0.2406138000566, 0.41002634575, 0.349359854193], 5:7)
testeval_snn(bs, 4.5, [0.0196263352654, 0.1775220289045, 0.80285163583], 5:7)
testeval_snn(bs, 5.0, [0.73141253905, 0.2318427624256, 0.0367446985247], 7:9)
testeval_snn(bs, 5.5, [0.291244696239, 0.414209710077, 0.294545593684], 7:9)
testeval_snn(bs, 6.0, [0.0355815798651, 0.228799037625, 0.73561938251], 7:9)
testeval_snn(bs, 2pi, [0, 0, 1], 7:9)
norms = sum(bs(collect(range(0, 2pi/99, 100)), coeffs) .^ 2, 2)
@test_approx_eq norms ones(100)

for j in 1:2
    for i in 0.5:0.5:6.0
        testeval_sdn(bs, i)
    end
    bs = deriv(bs)
end
@test_throws ArgumentError deriv(bs)
