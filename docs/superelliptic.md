# Superelliptic Curves

The `ascovers.superelliptic` package works with curves of the form

```text
C: y^m = f(x)
```

over the coefficient field of the univariate Sage polynomial `f`.  This part
of the package currently provides:

- `SuperellipticCurve`, with legacy alias `superelliptic`.
- `SuperellipticFunction`, with legacy alias `superelliptic_function`.
- `SuperellipticForm`, with legacy alias `superelliptic_form`.
- `SuperellipticDeRhamCocycle`, with legacy alias `superelliptic_cech`.
- Reduction helpers `reduction` and `reduction_form`.

The package is written as ordinary Python modules that use Sage objects, so the
examples below use `sage -python` compatible syntax.

## Imports

```python
from sage.all import GF, PolynomialRing
from ascovers import (
    SuperellipticCurve,
    SuperellipticDeRhamCocycle,
    SuperellipticFunction,
    SuperellipticForm,
    superelliptic,
    superelliptic_cech,
    superelliptic_function,
    superelliptic_form,
)
```

The lowercase names are compatibility aliases for the old Sage scripts.  New
code should prefer the capitalized class names.

## Constructing Curves

Create a curve from a univariate Sage polynomial and a positive integer
exponent:

```python
F = GF(7)
R = PolynomialRing(F, "x")
x = R.gen()

C = SuperellipticCurve(x**5 + x, 4)
```

This constructs the affine equation `y^4 = x^5 + x`.  The constructor
normalizes the defining polynomial into a polynomial ring over the coefficient
field and creates a two-variable coordinate ring with generators `x` and `y`.

The main constructor parameters are:

- `polynomial`: a univariate Sage polynomial over a field.
- `exponent`: the integer `m` in `y^m = f(x)`.
- `prec`: precision used for Laurent expansions at infinity.  The default is
  `100`.
- `compute_expansions`: whether to compute expansions at infinity during
  construction.  The default is `True`.
- `variable_name`: the name used for the normalized `x` variable.

In positive characteristic, `exponent` must be prime to the characteristic.
The defining polynomial must have positive degree.

## Basic Data

```python
C.polynomial
C.exponent
C.degree
C.base_ring
C.characteristic
C.nb_of_pts_at_infty
```

Useful curve methods include:

- `C.is_smooth()`: tests whether `f` is squarefree by checking the
  discriminant.
- `C.genus()`: returns
  `((deg(f) - 1) * (m - 1) - gcd(deg(f), m) + 1) / 2`.
- `C.infinity_parameters()`: returns the integer data used for the chart at
  infinity.
- `C.reciprocal_polynomial()`: returns `x^r f(1/x)`, where `r = deg(f)`.

The number of places at infinity is `delta = gcd(deg(f), m)`.  The current
expansion code works over the base field when the required roots at infinity
lie in the base field.  In particular, monic examples require a primitive
`delta`-th root of unity in the base field.

## Coordinate Functions and Forms

The curve provides convenient coordinate objects:

```python
C.x       # the function x
C.y       # the function y
C.one     # the constant function 1
C.dx      # the differential dx
```

These are instances of `SuperellipticFunction` and `SuperellipticForm`.
Arithmetic is reduced modulo the relation `y^m = f(x)`.

```python
F = GF(5)
R = PolynomialRing(F, "x")
x = R.gen()
C = SuperellipticCurve(x**3 + x + 1, 2)

C.y**2 == C.x**3 + C.x + C.one
```

You can also construct functions and forms from expressions in the coordinate
ring:

```python
x_coord, y_coord = C.coordinate_gens

g = C.function(x_coord + 2*y_coord + 1)
omega = C.form(y_coord)       # means y dx
```

The direct class constructors are equivalent:

```python
g = SuperellipticFunction(C, x_coord + 2*y_coord + 1)
omega = SuperellipticForm(C, y_coord)
```

## Reducing Expressions

The helper `reduction(C, expr)` reduces a rational expression in `x` and `y`
to a representative of degree less than `m` in `y`.

```python
from ascovers import reduction, reduction_form

reduced_function = reduction(C, y_coord**3 + x_coord*y_coord**2)
```

The helper `reduction_form(C, expr)` reduces a coefficient of a form to the
standard denominator form used by the differential-form class:

```text
sum_j h_j(x) dx / y^j
```

For example, on `y^2 = f(x)`, the form coefficient `y` is represented as
`f(x) / y` by this helper.  A `SuperellipticForm` then stores the coefficient
as a `SuperellipticFunction`, so its raw Sage representative may be any reduced
function equivalent modulo `y^m = f(x)`.

## Functions

`SuperellipticFunction` represents a rational function on `C`.

Common operations:

- Arithmetic: `+`, `-`, `*`, `/`, and powers.
- `g.jth_component(j)`: coefficient of `y^j` in the reduced representative.
- `g.diffn()`: returns the differential `dg` as a `SuperellipticForm`.
- `g.expansion_at_infty(place=0, prec=20)`: Laurent expansion at a place above
  infinity.
- `g.valuation(place=0)`: valuation at a place above infinity.
- `g.pth_root()`: in characteristic `p`, returns the `p`-th root when `g` is a
  `p`-th power, and raises `ValueError` otherwise.
- `g.numerator()` and `g.denominator()`: numerator and denominator as
  superelliptic functions.

Example:

```python
g = C.x + 2*C.y + C.one
dg = g.diffn()
```

Affine expansion at a point is available through `g.expansion(point, prec)`.
The current point interface is still minimal; pass an affine point as a pair
`(x0, y0)`.

The method `g.evaluate(point)` supports basic affine evaluation.  Evaluation at
infinity and cancellation through point objects will be completed when the
superelliptic point class is migrated.

## Differential Forms

`SuperellipticForm` represents a form `h(x, y) dx`.  Its coefficient is stored
as a `SuperellipticFunction` in `omega.form`; the underlying Sage
function-field expression is available as `omega.form.function`.

This is intentional: methods such as expansion and valuation reuse the
function implementation for the coefficient and then multiply by the expansion
of `dx`.

Common operations:

- Arithmetic: addition, subtraction, negation, and scalar/function
  multiplication.
- `omega.jth_component(j)`: the rational function `h_j(x)` when
  `omega = sum_j h_j(x) dx / y^j`.
- `omega.coordinates(basis=0)`: coordinates in a holomorphic differential
  basis.
- `omega.cartier()`: Cartier operator in positive characteristic.
- `omega.expansion_at_infty(place=0, prec=10)`: Laurent expansion at infinity.
- `omega.residue(place=0, prec=30)`: residue at infinity.
- `omega.valuation(place=0)`: valuation at infinity.
- `omega.is_regular_on_U0()` and `omega.is_regular_on_Uinfty()`: regularity
  checks on the two standard affine charts.
- `omega.serre_duality_pairing(g, prec=20)`: residue pairing with a function
  representing a class in `H^1(X, O_X)`.

Example:

```python
omega = C.y * C.dx
omega.form                    # SuperellipticFunction for y
omega.form.function           # raw Sage expression for the coefficient
omega.expansion_at_infty(prec=30)
```

## Cech-De Rham Cocycles

`SuperellipticDeRhamCocycle` represents a Cech-de Rham cocycle for the usual
two-open cover of a superelliptic curve.  It stores the triple

```text
(omega0, f, omega_infinity)
```

with the convention

```text
omega_infinity = omega0 - df.
```

The constructor takes the curve, the form `omega0`, and the transition
function `f`:

```python
eta = SuperellipticDeRhamCocycle(C, C.dx, C.x)
eta.omega0
eta.f
eta.omega8
eta.omega_infinity
```

The curve method `C.cech(omega0, f)` is equivalent.  The old name
`superelliptic_cech` is also available as a compatibility alias.

Common operations:

- Arithmetic: addition, subtraction, negation, and scalar multiplication.
- `eta.is_cocycle()`: returns `True` if `omega0` is regular on `U0` and
  `omega_infinity` is regular on the infinity chart; otherwise returns a short
  diagnostic string.
- `eta.is_valid()`: boolean wrapper around `eta.is_cocycle()`.
- `eta.verschiebung()` and `eta.frobenius()` in positive characteristic.
- `eta.coordinates(basis=0)`: coordinates in the standard de Rham basis for
  holomorphic, structure-sheaf, and mixed Cech representatives.

The standard de Rham basis is available from the curve:

```python
de_rham_basis = C.de_rham_basis()
for eta in de_rham_basis:
    eta.is_valid()
```

The helper `decomposition_g0_g8(g)` decomposes a function as `g0 - g8 + h`,
where `g0` is regular on the affine chart, `g8` is regular at infinity, and
`h` is a representative of the structure-sheaf cohomology part.

## Holomorphic Differentials

The method `C.holomorphic_differentials_basis()` returns the standard basis of
holomorphic differentials:

```python
basis = C.holomorphic_differentials_basis()
```

The basis consists of forms

```text
x^i dx / y^j
```

where the pairs `(i, j)` satisfy the holomorphicity inequality implemented by
the curve class.  To inspect the degree data:

```python
C.degrees_holomorphic_differentials()
```

Example:

```python
F = GF(7)
R = PolynomialRing(F, "x")
x = R.gen()
C = SuperellipticCurve(x**5 + x, 4)

C.degrees_holomorphic_differentials()
# {0: (0, 1), 1: (0, 2), 2: (1, 2), 3: (0, 3), 4: (1, 3), 5: (2, 3)}
```

Coordinates of a holomorphic form can be computed in this basis:

```python
basis = C.holomorphic_differentials_basis()
omega = 2*basis[0] + basis[1]
omega.coordinates(basis=basis)
```

## Laurent Expansions at Infinity

The constructor computes Laurent expansions of `x`, `y`, and `dx` at the
places above infinity unless `compute_expansions=False` is passed.

```python
C.x_series
C.y_series
C.dx_series
```

You normally use the object methods instead:

```python
C.x.expansion_at_infty(place=0, prec=50)
C.y.expansion_at_infty(place=0, prec=50)
(C.y * C.dx).expansion_at_infty(place=0, prec=50)
```

The valid places are indexed from `0` to `C.nb_of_pts_at_infty - 1`.

If a requested precision is higher than the stored precision, the curve
recomputes the expansions at that precision.

## Characteristic-p Operators

The current superelliptic implementation includes:

- `omega.cartier()` for differential forms.
- `C.cartier_matrix()` on holomorphic differentials.
- `C.a_number()`, computed as `genus - rank(Cartier matrix)`.
- `g.pth_root()` for functions that are `p`-th powers.
- `C.p_rank()` for hyperelliptic curves, using Sage's native
  `HyperellipticCurve.p_rank()`.

These methods require positive characteristic where appropriate.

Example:

```python
F = GF(67)
R = PolynomialRing(F, "x")
x = R.gen()
C = SuperellipticCurve(x**7 + x**3 + x, 2)

C.a_number()
```

## Compatibility with the Old Scripts

The old constructor names remain available:

```python
C = superelliptic(x**5 + x, 4)
g = superelliptic_function(C, C.coordinate_gens[0])
omega = superelliptic_form(C, C.coordinate_gens[1])
```

New code should use `SuperellipticCurve`, `SuperellipticFunction`, and
`SuperellipticForm`.  For de Rham cocycles, use
`SuperellipticDeRhamCocycle`; the old name `superelliptic_cech` remains
available.

## Current Migration Status

The curve, function, form, and Cech-de Rham cocycle classes are migrated.  Some methods on
`SuperellipticCurve` intentionally still depend on modules that have not been
migrated yet:

- General `SuperellipticDeRhamCocycle.coordinates()` needs the
  `decomposition_g0_g8` helper.  Coordinates of the standard de Rham basis and
  pure holomorphic or pure structure-sheaf cocycles are implemented.
- `C.final_type()` still needs cleanup of the Frobenius/Verschiebung matrix
  conventions.
- `C.riemann_roch_space()` and `C.riemann_roch_space_forms()` need the
  holomorphic-combination helpers from the cover code.
- `SuperellipticForm.inv_cartier()` needs the regular-form classes.

Those methods will become part of the stable documented API once the
corresponding legacy files have been migrated.

## Running Tests

From the repository root:

```bash
conda activate sage
DOT_SAGE=/tmp/ascovers-sage-dot PYTHONPATH=src sage -python -m unittest discover -s test -v
```
