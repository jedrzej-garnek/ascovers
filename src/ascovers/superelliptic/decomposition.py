'''Cech decomposition helpers for superelliptic curves.'''

from sage.all import FractionField, PolynomialRing


def decomposition_g0_g8(function, prec=50):
    '''Decompose a function as ``g0 - g8 + h`` for the two-chart cover.

    Here ``g0`` is regular on the affine chart, ``g8`` is regular at infinity,
    and ``h`` is a linear combination of the chosen basis of ``H^1(X, O_X)``.
    '''
    curve = function.curve
    coordinates = function.coordinates(prec=prec)
    nontrivial_part = 0 * curve.x
    for coefficient, basis_function in zip(coordinates, curve.cohomology_of_structure_sheaf_basis()):
        nontrivial_part += coefficient * basis_function

    function = function - nontrivial_part
    function_field, coordinate_ring, _x_coordinate, _y_coordinate = curve.fct_field
    rational_function = function_field(function.function)
    numerator = coordinate_ring(rational_function.numerator())
    denominator = coordinate_ring(rational_function.denominator())
    polynomial_part, remainder = numerator.quo_rem(denominator)

    denominator_function = curve.function(denominator)
    affine_part = curve.function(polynomial_part)
    infinity_part = 0 * curve.x

    for monomial in remainder.monomials():
        monomial_function = curve.function(monomial)
        quotient_regular_at_infinity = all(
            monomial_function.expansion_at_infty(place=place, prec=prec).valuation()
            >= denominator_function.expansion_at_infty(place=place, prec=prec).valuation()
            for place in range(int(curve.nb_of_pts_at_infty))
        )
        term = remainder.monomial_coefficient(monomial) * monomial_function / denominator_function
        if quotient_regular_at_infinity:
            infinity_part -= term
        else:
            affine_part += term

    return affine_part, infinity_part, nontrivial_part


def _coerce_univariate_to_curve_function(curve, rational_function):
    '''Evaluate a univariate rational function in the curve's x-coordinate.'''
    function_field, coordinate_ring, x_coordinate, _y_coordinate = curve.fct_field
    numerator = rational_function.numerator()
    denominator = rational_function.denominator()
    return curve.function(
        function_field(coordinate_ring(numerator(x_coordinate)))
        / function_field(coordinate_ring(denominator(x_coordinate)))
    )


def decomposition_omega0_omega8(form, prec=50):
    '''Decompose a residue-free form as ``omega0 - omega8`` on the two charts.'''
    curve = form.curve
    form = form.reduce()
    if sum(form.residue(place=place, prec=prec) for place in range(int(curve.nb_of_pts_at_infty))) != 0:
        raise ValueError(f"{form} has nonzero total residue")

    base_field = curve.base_ring
    polynomial_ring = PolynomialRing(base_field, "x")
    rational_field = FractionField(polynomial_ring)
    affine_function = 0 * curve.x
    infinity_function = 0 * curve.x

    for denominator_power in range(0, int(curve.exponent)):
        component = rational_field(form.jth_component(denominator_power))
        quotient, remainder = component.numerator().quo_rem(component.denominator())
        y_factor = curve.y ** (-denominator_power)
        affine_function += y_factor * curve.function(quotient(curve.coordinate_x))
        remainder_function = _coerce_univariate_to_curve_function(
            curve,
            rational_field(remainder) / rational_field(component.denominator()),
        )
        term = y_factor * remainder_function
        term_form = term * curve.dx
        if any(
            term_form.expansion_at_infty(place=place, prec=prec).valuation() < 0
            for place in range(int(curve.nb_of_pts_at_infty))
        ):
            raise ValueError(f"could not decompose {form}")
        infinity_function -= term

    affine_form = affine_function * curve.dx
    infinity_form = infinity_function * curve.dx
    if affine_form.is_regular_on_U0():
        return affine_form, infinity_form
    raise ValueError(f"could not decompose {form}; affine part is not regular")
