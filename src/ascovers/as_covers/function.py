'''Rational functions on Artin-Schreier type covers.'''

from sage.all import FractionField, Integer, LaurentSeriesRing, PolynomialRing, vector

from ascovers.as_covers.group import FiniteGroupElement
from ascovers.as_covers.reduction import artin_schreier_reduction


class ArtinSchreierFunction:
    '''A rational function on an Artin-Schreier type cover.'''

    def __init__(self, cover, expression):
        '''Create a function from a rational expression in ``x, y, z_i``.'''
        if isinstance(expression, ArtinSchreierFunction):
            if expression.curve is not cover:
                raise TypeError("cannot coerce a function from a different cover")
            expression = expression.function
        self.curve = cover
        function_field, _polynomial_ring, _x_coordinate, _y_coordinate, _z_coordinates = cover.fct_field
        self.function = function_field(expression)

    def __repr__(self):
        '''Return the textual representation of this function.'''
        return str(self.function)

    def __eq__(self, other):
        '''Return whether two functions agree on the cover.'''
        if not isinstance(other, ArtinSchreierFunction):
            try:
                other = ArtinSchreierFunction(self.curve, other)
            except (TypeError, ValueError):
                return False
        if other.curve is not self.curve:
            return False
        difference = self - other
        try:
            return difference.reduce().function == 0
        except Exception:
            series_ring = LaurentSeriesRing(self.curve.base_ring, "t", default_prec=self.curve.prec)
            expansion = series_ring(difference.expansion_at_infty())
            return expansion.valuation() >= 1

    def _coerce_function(self, other):
        '''Coerce another object to a function on the same cover.'''
        if isinstance(other, ArtinSchreierFunction):
            if other.curve is not self.curve:
                raise TypeError("cannot combine functions on different covers")
            return other
        return ArtinSchreierFunction(self.curve, other)

    def __add__(self, other):
        '''Add two functions or a function and a scalar expression.'''
        other = self._coerce_function(other)
        return ArtinSchreierFunction(self.curve, self.function + other.function)

    def __radd__(self, other):
        '''Add a scalar expression and this function.'''
        return self + other

    def __sub__(self, other):
        '''Subtract another function or scalar expression.'''
        other = self._coerce_function(other)
        return ArtinSchreierFunction(self.curve, self.function - other.function)

    def __rsub__(self, other):
        '''Subtract this function from a scalar expression.'''
        other = self._coerce_function(other)
        return ArtinSchreierFunction(self.curve, other.function - self.function)

    def __neg__(self):
        '''Return the additive inverse.'''
        return ArtinSchreierFunction(self.curve, -self.function)

    def __mul__(self, other):
        '''Multiply by a function, scalar expression, or AS differential form.'''
        from ascovers.as_covers.form import ArtinSchreierForm

        if isinstance(other, ArtinSchreierForm):
            if other.curve is not self.curve:
                raise TypeError("cannot multiply objects on different covers")
            return ArtinSchreierForm(self.curve, self * other.form)
        other = self._coerce_function(other)
        return ArtinSchreierFunction(self.curve, self.function * other.function)

    def __rmul__(self, other):
        '''Multiply this function by a scalar expression.'''
        return self * other

    def __truediv__(self, other):
        '''Divide by another function or scalar expression.'''
        other = self._coerce_function(other)
        return ArtinSchreierFunction(self.curve, self.function / other.function)

    def __rtruediv__(self, other):
        '''Divide a scalar expression by this function.'''
        other = self._coerce_function(other)
        return ArtinSchreierFunction(self.curve, other.function / self.function)

    def __pow__(self, exponent):
        '''Raise this function to an integer power.'''
        return ArtinSchreierFunction(self.curve, self.function ** exponent)

    def expansion_at_infty(self, place=0, prec=None):
        '''Return the Laurent expansion at a place above infinity.'''
        if prec is None:
            prec = self.curve.prec
        series_ring = LaurentSeriesRing(self.curve.base_ring, "t", default_prec=int(prec))
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.curve.fct_field
        place_key = int(place) if isinstance(place, (int, Integer)) else place
        substitutions = {
            x_coordinate: series_ring(self.curve.x_series[place_key]),
            y_coordinate: series_ring(self.curve.y_series[place_key]),
        }
        substitutions.update(
            {z_coordinates[index]: series_ring(self.curve.z_series[place_key][index]) for index in range(self.curve.height)}
        )
        return series_ring(function_field(self.function).substitute(substitutions))

    def expansion(self, point=0, prec=None):
        '''Return a Laurent expansion at a currently supported branch point.'''
        return self.expansion_at_infty(place=point, prec=prec)

    def group_action(self, element):
        '''Return the image of this function under a group element.'''
        cover = self.curve
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = cover.fct_field
        template_field, template_z, template_f, template_x, template_y = cover.cover_template.fct_field

        if isinstance(element, FiniteGroupElement):
            element = element.as_tuple

        if element == cover.group.one:
            return self

        if element in cover.group.gens:
            generator_index = cover.group.gens.index(element)
            action = cover.cover_template.gp_action[generator_index]
            template_substitutions = {
                template_z[index]: function_field(z_coordinates[index])
                for index in range(cover.height)
            }
            template_substitutions.update(
                {
                    template_f[index]: function_field(cover.functions[index].function)
                    for index in range(cover.height)
                }
            )
            template_substitutions.update({template_x: x_coordinate, template_y: y_coordinate})
            substitutions = {
                x_coordinate: function_field(action[-2].substitute(template_substitutions)),
                y_coordinate: function_field(action[-1].substitute(template_substitutions)),
            }
            substitutions.update(
                {
                    z_coordinates[index]: function_field(action[index].substitute(template_substitutions))
                    for index in range(cover.height)
                }
            )
            return ArtinSchreierFunction(cover, self.function.substitute(substitutions))

        result = self
        if isinstance(element, (list, tuple)):
            exponents = element
        else:
            exponents = [element] * len(cover.group.gens)
        for generator, exponent in zip(cover.group.gens, exponents):
            for _step in range(int(exponent)):
                result = result.group_action(generator)
        return result

    def reduce(self):
        '''Return this function reduced modulo the cover equations.'''
        return ArtinSchreierFunction(self.curve, artin_schreier_reduction(self.curve, self.function))

    def trace(self, super=True, subgp=0):
        '''Return the trace of this function over the requested subgroup.

        When ``super`` is true and the full group trace is used, the result is
        returned as a function on the quotient superelliptic curve.
        '''
        cover = self.curve
        if isinstance(subgp, (int, Integer)):
            subgroup = cover.group.elts
        else:
            subgroup = subgp
            super = False

        result = 0 * cover.x
        for element in subgroup:
            result += self.group_action(element)
        reduced = artin_schreier_reduction(cover, result.function)

        if super:
            from ascovers.superelliptic.function import SuperellipticFunction

            quotient_ring = PolynomialRing(cover.base_ring, ("x", "y"))
            quotient_field = FractionField(quotient_ring)
            return SuperellipticFunction(cover.quotient, quotient_field(reduced))
        return ArtinSchreierFunction(cover, reduced)

    def coordinates(self, prec=100, basis=0):
        '''Return coordinates in H^1(X, O_X), dual to holomorphic forms.'''
        cover = self.curve
        if basis == 0:
            basis = [
                cover.holomorphic_differentials_basis(),
                cover.cohomology_of_structure_sheaf_basis(),
            ]
        holomorphic_forms, cohomology_basis = basis
        pairing_matrix = [
            [omega.serre_duality_pairing(function) for omega in holomorphic_forms]
            for function in cohomology_basis
        ]
        target = [omega.serre_duality_pairing(self) for omega in holomorphic_forms]
        vector_space = (cover.base_ring ** int(cover.genus())).span_of_basis(
            [vector(cover.base_ring, row) for row in pairing_matrix]
        )
        return vector_space.coordinates(vector(cover.base_ring, target))

    def diffn(self):
        '''Return the differential of this function as a multiple of dx.'''
        cover = self.curve
        quotient = cover.quotient
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = cover.fct_field
        function = function_field(self.function)
        dy_coefficient = function_field(quotient.y.diffn().form.function)
        result = function.derivative(x_coordinate)
        result += function.derivative(y_coordinate) * dy_coefficient
        for index in range(cover.height):
            result += function.derivative(z_coordinates[index]) * cover.dz[index]

        from ascovers.as_covers.form import ArtinSchreierForm

        return ArtinSchreierForm(cover, result)

    def valuation(self, place=0):
        '''Return the valuation at a place above infinity.'''
        series_ring = LaurentSeriesRing(self.curve.base_ring, "t", default_prec=self.curve.prec)
        return series_ring(self.expansion_at_infty(place=place)).valuation()

    def numerator(self):
        '''Return the numerator as a function on the same cover.'''
        return ArtinSchreierFunction(self.curve, self.function.numerator())

    def denominator(self):
        '''Return the denominator as a function on the same cover.'''
        return ArtinSchreierFunction(self.curve, self.function.denominator())


as_function = ArtinSchreierFunction
