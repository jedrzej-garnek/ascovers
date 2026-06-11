'''Utility functions for computations with superelliptic curves.'''

from sage.all import FractionField, LaurentSeriesRing, PolynomialRing, ZZ


def coefficient(polynomial, degree):
    '''Return the coefficient of a univariate polynomial at the requested degree.'''
    if degree < 0:
        return polynomial.parent().zero()
    return polynomial[degree]


def cut_polynomial(polynomial, index):
    '''Return the polynomial formed from terms of degree strictly larger than index.'''
    polynomial_ring = polynomial.parent()
    variable = polynomial_ring.gen()
    result = polynomial_ring.zero()
    for degree in range(index + 1, polynomial.degree() + 1):
        result += coefficient(polynomial, degree) * variable ** (degree - index - 1)
    return result


def cut(polynomial, index):
    '''Compatibility alias for the legacy helper name.'''
    return cut_polynomial(polynomial, index)


def degree_of_rational_function(function, base_field):
    '''Return degree(numerator) minus degree(denominator) for a rational function.'''
    polynomial_ring = PolynomialRing(base_field, "x")
    rational_function_field = FractionField(polynomial_ring)
    rational_function = rational_function_field(function)
    return rational_function.numerator().degree() - rational_function.denominator().degree()


def coefficient_of_rational_function(function, base_field):
    '''Return the quotient of leading coefficients of numerator and denominator.'''
    polynomial_ring = PolynomialRing(base_field, "x")
    rational_function_field = FractionField(polynomial_ring)
    rational_function = rational_function_field(function)
    if rational_function == polynomial_ring.zero():
        return base_field.zero()
    numerator = rational_function.numerator()
    denominator = rational_function.denominator()
    return numerator.leading_coefficient() / denominator.leading_coefficient()


def polynomial_part(characteristic, polynomial):
    '''Return the Cartier polynomial part with exponents congruent to p - 1 modulo p.'''
    base_field = polynomial.parent().base_ring()
    polynomial_ring = PolynomialRing(base_field, "x")
    variable = polynomial_ring.gen()
    polynomial = polynomial_ring(polynomial)
    result = polynomial_ring.zero()
    for degree in range(polynomial.degree() + 1):
        if degree % characteristic == characteristic - 1:
            result += base_field(polynomial[degree]) * variable ** ((degree - characteristic + 1) // characteristic)
    return result


def root_of_unity(field, order):
    '''Return a primitive root of unity of the requested order in a field.'''
    order = ZZ(order)
    if order <= 0:
        raise ValueError("the order of a root of unity must be positive")
    if order == 1:
        return field.one()

    polynomial_ring = PolynomialRing(field, "z")
    variable = polynomial_ring.gen()
    for root, _multiplicity in (variable ** order - 1).roots():
        proper_powers = [root ** divisor for divisor in order.divisors() if divisor != order]
        if field.one() not in proper_powers:
            return root

    if hasattr(field, "base_ring"):
        base_ring = field.base_ring()
        if base_ring is not field:
            return field(root_of_unity(base_ring, order))

    raise ValueError(f"the field does not contain a primitive {order}-th root of unity")


def roots_in_base_field(field, exponent, value):
    '''Return all roots in the field of z^exponent - value.'''
    exponent = ZZ(exponent)
    if exponent <= 0:
        raise ValueError("the exponent must be positive")
    polynomial_ring = PolynomialRing(field, "z")
    variable = polynomial_ring.gen()
    roots = [root for root, _multiplicity in (variable ** exponent - field(value)).roots()]
    if len(roots) != exponent:
        raise ValueError(
            f"z^{exponent} - ({value}) does not split into {exponent} roots over the base field"
        )
    return roots


def naive_hensel(polynomial, field, start=1, prec=10):
    '''Lift a simple root of a polynomial over F((t)) by Newton-Hensel iteration.'''
    precision = int(prec)
    series_ring = LaurentSeriesRing(field, "t", default_prec=precision)
    series_field = FractionField(series_ring)
    polynomial_ring = PolynomialRing(series_field, "W")
    polynomial_field = FractionField(polynomial_ring)

    hensel_polynomial = polynomial_ring(polynomial_field(polynomial).numerator())
    derivative_at_start = hensel_polynomial.derivative()(series_field(start))
    if derivative_at_start == 0:
        raise ZeroDivisionError("the starting value is not a simple root modulo t")

    approximation = series_ring(start)
    for _step in range(1, precision):
        correction = hensel_polynomial(approximation) / derivative_at_start
        approximation = series_ring(approximation - correction).add_bigoh(precision)
    return approximation


def nth_root_series(series, exponent, prec=10):
    '''Return an n-th root of a Laurent series when the chosen root exists.'''
    exponent = ZZ(exponent)
    field = series.parent().base_ring()
    precision = int(prec)
    series_ring = LaurentSeriesRing(field, "t", default_prec=precision)
    variable = series_ring.gen()
    polynomial_ring = PolynomialRing(series_ring, "W")
    polynomial_variable = polynomial_ring.gen()

    valuation = series.valuation()
    if valuation % exponent != 0:
        raise ValueError("the valuation of the series is not divisible by the exponent")

    normalized_series = series_ring(series * variable ** (-valuation))
    constant_term = normalized_series[0]
    lifted_root = naive_hensel(
        polynomial_variable ** exponent - normalized_series,
        field,
        start=constant_term.nth_root(exponent),
        prec=precision,
    )
    return variable ** (valuation // exponent) * lifted_root


def reverse_power_series(power_series, prec=10):
    '''Return the compositional inverse of a power series with nonzero linear term.'''
    field = power_series.parent().base_ring()
    precision = int(prec)
    series_ring = LaurentSeriesRing(field, "t", default_prec=precision)
    series_field = FractionField(series_ring)
    variable = series_ring.gen()
    power_series = series_field(power_series)

    linear_coefficient = power_series.list()[0]
    inverse_series = (1 / linear_coefficient) * variable
    for degree in range(2, precision + 1):
        residual = power_series(inverse_series) - variable
        coefficient = 0 if residual.valuation() > degree else residual.list()[0]
        inverse_series -= coefficient / linear_coefficient * variable ** degree
    return inverse_series


def new_reverse(power_series, prec=10):
    '''Compatibility alias for the legacy compositional-inverse helper.'''
    return reverse_power_series(power_series, prec=prec)
