'''Reduction of rational functions on Artin-Schreier type covers.'''

from sage.all import prod


def artin_schreier_reduction(cover, expression):
    '''Reduce a rational expression using the cover equations.

    The result is represented with powers of each Artin-Schreier generator below
    ``p`` and with the superelliptic relation on the quotient applied to powers of
    ``y``.  Denominators involving AS generators are first cleared by multiplying
    by their group norm.
    '''
    from ascovers.as_covers.function import ArtinSchreierFunction

    function_field, polynomial_ring, x_coordinate, y_coordinate, z_coordinates = cover.fct_field
    expression = function_field(expression)
    numerator = expression.numerator()
    denominator = expression.denominator()

    if denominator != 1:
        denominator_function = ArtinSchreierFunction(cover, denominator)
        denominator_norm = cover.one
        for group_element in cover.group.elts:
            if group_element != cover.group.one:
                denominator_norm *= denominator_function.group_action(group_element)
        reduced_numerator = artin_schreier_reduction(
            cover, polynomial_ring(numerator * denominator_norm.function)
        )
        reduced_denominator = artin_schreier_reduction(
            cover, polynomial_ring(denominator * denominator_norm.function)
        )
        return function_field(reduced_numerator) / function_field(reduced_denominator)

    quotient_polynomial = polynomial_ring(cover.quotient.polynomial(x_coordinate))
    quotient_exponent = int(cover.quotient.exponent)
    characteristic = int(cover.characteristic)
    numerator = polynomial_ring(numerator)
    result = function_field.zero()
    changed = False

    for monomial in numerator.monomials():
        z_degrees = [monomial.degree(z_coordinate) for z_coordinate in z_coordinates]
        z_quotients = [degree // characteristic for degree in z_degrees]
        z_remainders = [degree % characteristic for degree in z_degrees]
        y_degree = monomial.degree(y_coordinate)
        if any(z_quotients) or y_degree >= quotient_exponent:
            changed = True

        coefficient = numerator.monomial_coefficient(monomial)
        reduced_monomial = (
            coefficient
            * x_coordinate ** monomial.degree(x_coordinate)
            * y_coordinate ** (y_degree % quotient_exponent)
            * quotient_polynomial ** (y_degree // quotient_exponent)
            * prod(
                z_coordinates[index] ** z_remainders[index]
                for index in range(len(z_coordinates))
            )
            * prod(
                (z_coordinates[index] + cover.rhs[index]) ** z_quotients[index]
                for index in range(len(z_coordinates))
            )
        )
        result += function_field(reduced_monomial)

    if not changed:
        return function_field(result)
    return artin_schreier_reduction(cover, result)


as_reduction = artin_schreier_reduction
