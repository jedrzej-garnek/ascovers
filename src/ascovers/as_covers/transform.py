'''Local transformations for Artin-Schreier equations.'''

from sage.all import FractionField, Infinity, LaurentSeriesRing

from ascovers.superelliptic.utils import nth_root_series, reverse_power_series


def _coefficient_field(power_series):
    '''Return the coefficient field of a Laurent series or its fraction field.'''
    parent = power_series.parent()
    if hasattr(parent, "base_ring"):
        base = parent.base_ring()
    else:
        base = parent.base()
    if hasattr(base, "base_ring") and "Laurent Series" in str(base):
        return base.base_ring()
    return base


def artin_schreier_transform(power_series, prec=10):
    '''Normalize a local Artin-Schreier equation and compute the new parameter.

    Given a Laurent series ``f(t)``, this returns data for the local equation
    ``z^p - z = f(t)`` after repeatedly removing terms whose negative valuation
    is divisible by ``p``.  The returned tuple is ``(jump, correction, t_old, z)``,
    where ``t_old`` expresses the old parameter in the new local parameter and
    ``z`` is the corresponding Laurent expansion of the Artin-Schreier generator.
    '''
    coefficient_field = _coefficient_field(power_series)
    characteristic = coefficient_field.characteristic()
    series_ring = LaurentSeriesRing(coefficient_field, "t", default_prec=int(prec))
    series_field = FractionField(series_ring)
    parameter = series_ring.gen()
    normalized_rhs = series_field(power_series)
    correction = series_ring.zero()

    if normalized_rhs.valuation() == Infinity:
        raise ValueError("precision is too low to normalize the Artin-Schreier equation")

    if normalized_rhs.valuation() >= 0:
        inverse = reverse_power_series(parameter ** characteristic - parameter, prec=prec)
        return 0, correction, parameter, inverse(normalized_rhs)

    while normalized_rhs.valuation() < 0 and normalized_rhs.valuation() % characteristic == 0:
        pole_order = int(-normalized_rhs.valuation() // characteristic)
        leading_coefficient = normalized_rhs.list()[0]
        root = leading_coefficient.nth_root(characteristic)
        correction += root * parameter ** (-pole_order)
        normalized_rhs -= leading_coefficient * parameter ** (-characteristic * pole_order)
        normalized_rhs += root * parameter ** (-pole_order)

    jump = max(int(-normalized_rhs.valuation()), 0)
    if jump == 0:
        inverse = reverse_power_series(parameter ** characteristic - parameter, prec=prec)
        generator_series = inverse(normalized_rhs) + correction
        return 0, correction, parameter, generator_series

    try:
        scaled_parameter = nth_root_series(normalized_rhs ** (-1), jump, prec=prec)
        denominator = nth_root_series(1 - parameter ** ((characteristic - 1) * jump), jump, prec=prec)
    except Exception as error:
        raise ValueError(f"could not extract the {jump}-th root needed for the local transform") from error

    reverse_scaled_parameter = reverse_power_series(scaled_parameter, prec=prec)
    old_parameter = reverse_scaled_parameter(parameter ** characteristic / denominator)
    generator_series = parameter ** (-jump) + series_ring(correction)(old_parameter)
    return jump, correction, old_parameter, generator_series


as_transform = artin_schreier_transform
