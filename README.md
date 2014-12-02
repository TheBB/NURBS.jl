# NURBS

[![Build Status](https://travis-ci.org/TheBB/NURBS.jl.svg?branch=master)](https://travis-ci.org/TheBB/NURBS.jl)

**This package does not support the latest Julia release version (0.3).  You need to use the 0.4
  development version.**

_NURBS.jl_ is a Julia package for manipulating NURBS objects (curves, surfaces, volumes, etc.)  This
includes B-Splines, as a special case.  It involves no external dependencies.

This is written mostly as a way to teach myself Julia.  Other packages provide similar functionality:

* [BSplines](https://github.com/gusl/BSplines.jl) by gusl. Failing tests as of the time of writing.
* [Dierckx](https://github.com/kbarbary/Dierckx.jl) by kbarbary. Has external dependencies.
* [Grid](https://github.com/timholy/Grid.jl) by timholy. Only deals with uniformly spaced splines.
* [Interpolations](https://github.com/tlycken/Interpolations.jl) by tlycken. A continuation of _Grid_.

_NURBS.jl_ intends to provide a common framework objects of both spline and NURBS types, with no
external dependencies.

## Status

The package is currently very early in development.  See below for what is implemented.

## Usage

### Creating basis objects

Import `NURBS.Bases`.  The two bases currently implemented are `BSplineBasis` and `NURBSBasis`, both
subtypes of `Basis1D`.

Instantiate a `BSplineBasis` with a knot vector and a given order:

```julia
basis = BSplineBasis([0, 1, 3, 4], 4)
```

Note that the order is one higher than the polynomial degree of the basis functions.  An order 4
basis is cubic.  The constructor will create repeated knots at the endpoints automatically.  Pass
`extend=false` as a keyword argument if you don't want this.

You can also create a uniform basis more easily.

```julia
basis = BSplineBasis(lower, upper, nelements, order)
    # Equivalent to
    # basis = BSplineBasis(linspace(lower, upper, nelements+1), order)
```

Note that, unlike `linspace`, this constructor expects the number of _elements_ (intervals), not the
number of _points_.

A NURBS basis is instantiated from a B-Spline basis and a vector of weights, one for each basis
function.

```julia
nbasis = NURBSBasis(bbasis, [1.0, 2.0, ...])
```

### Using basis objects

* `length` and `size` give the basis dimensions, much like how they work with arrays.  Note that
  since only one-dimensional bases are implemented currently, you will always have `size(basis) ==
  (length(basis),)`.
* `domain` gives the interval on which a basis is supported as an `Interval` object which supports
  the operators `in`, `==`, ⊆, ⊈ and ⊊, but _not_ iteration.
* `degree` gives the polynomial degree of basis functions.
* `nderivs` gives the number of supported derivatives for a basis (this should be equal to the
  degree for B-Splines).
* `deriv` returns a _new_ basis object whose functions are the derivatives of those in the original,
  provided the original object supports differentiation (`nderivs()` > 0).  Differentiate by
  multiple orders with `deriv(basis, n)`.
* You can also index a basis object to obtain a basis function.  Basis functions support the
  `domain`, `degree` and `deriv` operations.

### Evaluating basis functions

_NURBS.jl_ basis and basis function objects make use of Julia's call overloading, which is present
in 0.4 but not in earlier version.  For this reason, _NURBS.jl_ will not work on Julia 0.3.

Evaluation at a single point (e.g. `basis(1)`) will return a tuple `(vals, idxs)`, where `vals` is a
vector of basis function values at that point, and `idxs` is a range denoting which basis functions
are supported there.  Only the supported basis functions are actually evaluated, although this
doesn't mean they cannot be zero there!

You can evaluate a single basis function at a single point, which will return a `Float64`.

Evaluation at multiple points will return, for a basis, an array of tuples, or for a basis function,
an array of floats.

You can also evaluate a basis with coefficients.  Coefficients are stored either in a vector of size
`length(basis)` or a matrix of size `length(basis)` × _d_, where _d_ is the dimension of the
surrounding space.  Calls to `basis(pt, coeffs)` will return a float or a 1×_d_-dimensional matrix,
while calls to `basis(pts, coeffs)` will return a vector or a matrix.

(This may change when I decide it is more convenient to do something else.)

### Example

The unit circle in NURBS form.  See the
[Wikipedia article](http://en.wikipedia.org/wiki/Non-uniform_rational_B-spline#Example:_a_circle).

```julia
julia> using NURBS.Bases

julia> knots = [0, pi/2, pi/2, pi, pi, 3pi/2, 3pi/2, 2pi];

julia> bsplines = BSplineBasis(knots, 3);

julia> t = sqrt(2) / 2;

julia> weights = [1, t, 1, t, 1, t, 1, t, 1];

julia> nurbs = NURBSBasis(bsplines, weights);

julia> coeffs = [1 0; 1 1; 0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1; 1 0];

julia> for p in linspace(0, 2pi, 9) @show nurbs(p, coeffs) end
nurbs(p,coeffs) = [1.0 0.0]
nurbs(p,coeffs) = [0.7071067811865476 0.7071067811865476]
nurbs(p,coeffs) = [0.0 1.0]
nurbs(p,coeffs) = [-0.7071067811865476 0.7071067811865476]
nurbs(p,coeffs) = [-1.0 0.0]
nurbs(p,coeffs) = [-0.7071067811865476 -0.7071067811865476]
nurbs(p,coeffs) = [0.0 -1.0]
nurbs(p,coeffs) = [0.7071067811865476 -0.7071067811865476]
nurbs(p,coeffs) = [1.0 0.0]

julia> points = nurbs(linspace(0, 2pi, 1000), coeffs);

julia> (min, max) = extrema(sum(points .^ 2, 2))
(0.9999999999999997,1.0000000000000004)
```
