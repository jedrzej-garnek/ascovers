'''Differential forms on Artin-Schreier type covers.'''

from sage.all import FractionField, Integer, LCM, LaurentSeriesRing, PolynomialRing, prod

from ascovers.as_covers.reduction import artin_schreier_reduction
from ascovers.linear_algebra import linear_representation_polynomials


class ArtinSchreierForm:
    '''A differential form on an AS cover, represented as ``h dx``.'''

    def __init__(self, cover, expression):
        '''Create the form represented by ``expression * dx``.'''
        from ascovers.as_covers.function import ArtinSchreierFunction

        if isinstance(expression, ArtinSchreierForm):
            if expression.curve is not cover:
                raise TypeError("cannot coerce a form from a different cover")
            coefficient = expression.form
        elif isinstance(expression, ArtinSchreierFunction):
            if expression.curve is not cover:
                raise TypeError("cannot coerce a function from a different cover")
            coefficient = expression
        else:
            coefficient = ArtinSchreierFunction(cover, expression)
        self.curve = cover
        self.form = coefficient

    def __repr__(self):
        '''Return the textual representation of this form.'''
        if len(str(self.form)) == 1:
            return str(self.form) + " dx"
        return "(" + str(self.form) + ") dx"

    def __eq__(self, other):
        '''Return whether two forms have the same coefficient on the cover.'''
        return isinstance(other, ArtinSchreierForm) and self.curve is other.curve and self.form == other.form

    def __add__(self, other):
        '''Add two differential forms.'''
        if other.curve is not self.curve:
            raise TypeError("cannot add forms on different covers")
        return ArtinSchreierForm(self.curve, self.form + other.form)

    def __sub__(self, other):
        '''Subtract another differential form.'''
        if other.curve is not self.curve:
            raise TypeError("cannot subtract forms on different covers")
        return ArtinSchreierForm(self.curve, self.form - other.form)

    def __neg__(self):
        '''Return the additive inverse.'''
        return ArtinSchreierForm(self.curve, -self.form)

    def __rmul__(self, constant):
        '''Multiply this form by a scalar or AS function.'''
        from ascovers.as_covers.function import ArtinSchreierFunction

        if isinstance(constant, ArtinSchreierFunction):
            if constant.curve is not self.curve:
                raise TypeError("cannot multiply objects on different covers")
            return ArtinSchreierForm(self.curve, constant * self.form)
        return ArtinSchreierForm(self.curve, constant * self.form)

    def expansion_at_infty(self, place=0, prec=None):
        '''Return the Laurent expansion at a place above infinity.'''
        if prec is None:
            prec = self.curve.prec
        series_ring = LaurentSeriesRing(self.curve.base_ring, "t", default_prec=int(prec))
        place_key = int(place) if isinstance(place, (int, Integer)) else place
        coefficient = self.form.expansion_at_infty(place=place_key, prec=prec)
        return series_ring(coefficient) * series_ring(self.curve.dx_series[place_key])

    def expansion(self, point=0, prec=None):
        '''Return a Laurent expansion at a currently supported branch point.'''
        return self.expansion_at_infty(place=point, prec=prec)

    def reduce(self):
        '''Return this form with its coefficient reduced modulo the cover equations.'''
        return ArtinSchreierForm(
            self.curve,
            artin_schreier_reduction(self.curve, self.form.function),
        )

    def group_action(self, element):
        '''Return the image of this form under a group element.'''
        dx_action = self.curve.x.group_action(element).diffn()
        coefficient_action = self.form.group_action(element)
        return coefficient_action * dx_action

    def coordinates(self, basis=0):
        '''Return coordinates of this holomorphic form in a given basis.'''
        reduced_form = self.reduce()
        cover = self.curve
        if basis == 0:
            basis = cover.holomorphic_differentials_basis()
        _function_field, polynomial_ring, _x_coordinate, _y_coordinate, _z_coordinates = cover.fct_field
        common_denominator = LCM([omega.form.function.denominator() for omega in basis])
        basis_without_denominators = [common_denominator * omega.form.function for omega in basis]
        form_without_denominator = common_denominator * reduced_form.form.function
        return linear_representation_polynomials(
            polynomial_ring(form_without_denominator),
            [polynomial_ring(omega) for omega in basis_without_denominators],
        )

    def trace(self, super=True):
        '''Return the trace of this form over the full Galois group.'''
        cover = self.curve
        result = 0 * cover.dx
        for element in cover.group.elts:
            result += self.group_action(element)
        reduced = artin_schreier_reduction(cover, result.form.function)
        if super:
            from ascovers.superelliptic.form import SuperellipticForm

            quotient_ring = PolynomialRing(cover.base_ring, ("x", "y"))
            quotient_field = FractionField(quotient_ring)
            return SuperellipticForm(cover.quotient, quotient_field(reduced))
        return ArtinSchreierForm(cover, reduced)

    def residue(self, place=0, prec=None):
        '''Return the residue at a place above infinity.'''
        return self.expansion_at_infty(place=place, prec=prec)[-1]

    def residues(self, place=0, prec=None):
        '''Return residues of all group translates at a place above infinity.'''
        return [self.group_action(element).residue(place=place, prec=prec) for element in self.curve.group.elts]

    def valuation(self, place=0):
        '''Return the valuation at a place above infinity.'''
        return self.expansion_at_infty(place=place).valuation()

    def serre_duality_pairing(self, function):
        '''Pair this form with a class in H^1(X, O_X) by residues.'''
        return sum((function * self).residue(place=place) for place in range(int(self.curve.nb_of_pts_at_infty)))

    def cartier(self):
        '''Compute the Cartier operator on this form when it is in quotient form.'''
        cover = self.curve
        base_field = cover.base_ring
        characteristic = base_field.characteristic()
        quotient = cover.quotient
        function_field, polynomial_ring, x_coordinate, y_coordinate, z_coordinates = cover.fct_field
        coefficient = function_field(self.form.function)
        quotient_ring = PolynomialRing(base_field, ("x", "y"))
        quotient_field = FractionField(quotient_ring)

        numerator = polynomial_ring(coefficient.numerator())
        denominator = polynomial_ring(coefficient.denominator())
        result = function_field.zero()
        if denominator not in quotient_ring:
            raise ValueError(
                "present the form as a sum z^i omega_i with omega_i forms on the quotient curve"
            )

        substitutions = {
            x_coordinate: x_coordinate,
            y_coordinate: y_coordinate,
        }
        substitutions.update(
            {
                z_coordinates[index]: z_coordinates[index] ** characteristic
                - function_field(cover.functions[index].function)
                for index in range(cover.height)
            }
        )
        expanded_numerator = function_field(numerator.substitute(substitutions))
        auxiliary_denominator = polynomial_ring(expanded_numerator.denominator())
        expanded_numerator = polynomial_ring(expanded_numerator * auxiliary_denominator ** characteristic)

        from ascovers.superelliptic.form import SuperellipticForm

        for monomial in expanded_numerator.monomials():
            z_degrees = [monomial.degree(z_coordinates[index]) for index in range(cover.height)]
            z_part = prod(
                z_coordinates[index] ** z_degrees[index]
                for index in range(cover.height)
            )
            quotient_part = monomial / z_part
            z_part_without_powers = prod(
                z_coordinates[index] ** (z_degrees[index] // characteristic)
                for index in range(cover.height)
            )
            quotient_form = SuperellipticForm(quotient, quotient_field(quotient_part / denominator))
            quotient_cartier = quotient_form.cartier()
            result += (
                z_part_without_powers
                * expanded_numerator.monomial_coefficient(monomial)
                * function_field(quotient_cartier.form.function)
                / auxiliary_denominator
            )
        return ArtinSchreierForm(cover, result)

    def is_regular_on_U0(self):
        '''Return whether this form is regular on the affine chart.'''
        cover = self.curve
        _function_field, polynomial_ring, _x_coordinate, y_coordinate, _z_coordinates = cover.fct_field
        return y_coordinate ** (cover.quotient.exponent - 1) * self.form.function in polynomial_ring


def are_forms_linearly_dependent(forms):
    '''Return whether the given AS forms are linearly dependent over the base field.'''
    from sage.rings.polynomial.toy_variety import is_linearly_dependent

    cover = forms[0].curve
    _function_field, polynomial_ring, _x_coordinate, _y_coordinate, _z_coordinates = cover.fct_field
    common_denominator = prod(form.form.function.denominator() for form in forms)
    return is_linearly_dependent(
        [polynomial_ring(common_denominator * form.form.function) for form in forms]
    )


def only_log_forms(cover):
    '''Return forms with at most logarithmic poles not generated by holomorphic forms.'''
    holomorphic_forms = cover.at_most_poles_forms(0)
    logarithmic_forms = cover.at_most_poles_forms(1)
    result = []
    for form in logarithmic_forms:
        if not are_forms_linearly_dependent(holomorphic_forms + result + [form]):
            result.append(form)
    return result


as_form = ArtinSchreierForm
