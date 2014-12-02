using NURBS.Bases

bs = BSplineBasis(0, 1, 1, 4)
pts = linspace(0, 1, 100000)

bs(pts)
@time bs(pts)
