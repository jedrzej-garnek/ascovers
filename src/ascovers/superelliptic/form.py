'''Differential forms on superelliptic curves.'''

from sage.all import (
    FractionField,
    Integers,
    LCM,
    LaurentSeriesRing,
    PolynomialRing,
    ZZ,
    floor,
)

from ascovers.linear_algebra import linear_representation_polynomials
from ascovers.superelliptic.curve import _expression_to_polynomial_in_y, reduction, reduction_form
from ascovers.superelliptic.exceptions import PendingMigrationError
from ascovers.superelliptic.utils import degree_of_rational_function, polynomial_part


def _coerce_univariate_rational(rational_function, target_field):
    '''Coerce a univariate rational function into another rational function field.'''
    target_ring = target_field.ring()
    target_variable = target_ring.gen()
    numerator = rational_function.numerator()
    denominator = rational_function.denominator()
    return target_field(target_ring(numerator(target_variable))) / target_field(
        target_ring(denominator(target_variable))
    )


class SuperellipticForm:
    '''A differential form on a superelliptic curve, represented as h(x, y) dx.'''

    def __init__(self, curve, expression):
        '''Create the differential form represented by expression times dx.'''
        self.curve = curve
        self.form = curve.function_field(reduction_form(curve, expression))

    def __eq__(self, other):
        '''Return whether two forms have the same reduced coefficient.'''
        if not isinstance(other, SuperellipticForm):
            return False
        return self.curve is other.curve and self.reduce().form == other.reduce().form

    def __add__(self, other):
        '''Add two differential forms.'''
        if other.curve is not self.curve:
            raise TypeError("cannot add forms on different curves")
        return SuperellipticForm(self.curve, self.form + other.form)

    def __sub__(self, other):
        '''Subtract another differential form.'''
        if other.curve is not self.curve:
            raise TypeError("cannot subtract forms on different curves")
        return SuperellipticForm(self.curve, self.form - other.form)

    def __neg__(self):
        '''Return the additive inverse of this form.'''
        return SuperellipticForm(self.curve, -self.form)

    def __repr__(self):
        '''Return the textual representation of this form.'''
        if len(str(self.form)) == 1:
            return str(self.form) + " dx"
        return "(" + str(self.form) + ") dx"

    def __rmul__(self, other):
        '''Multiply this form by a scalar or superelliptic function.'''
        from ascovers.superelliptic.function import SuperellipticFunction

        if isinstance(other, SuperellipticFunction):
            if other.curve is not self.curve:
                raise TypeError("cannot multiply objects on different curves")
            return SuperellipticForm(self.curve, other.function * self.form)
        return SuperellipticForm(self.curve, other * self.form)

    def __truediv__(self, other):
        '''Return the quotient of coefficients of two differential forms.'''
        from ascovers.superelliptic.function import SuperellipticFunction

        if other.curve is not self.curve:
            raise TypeError("cannot divide forms on different curves")
        return SuperellipticFunction(self.curve, self.form / other.form)

    def cartier(self):
        '''Compute the Cartier operator on this differential form.'''
        characteristic = self.curve.characteristic
        if characteristic == 0:
            raise ValueError("Cartier operator is only implemented in positive characteristic")

        base_field = self.curve.base_ring
        polynomial_ring = PolynomialRing(base_field, "x")
        x_variable = polynomial_ring.gen()
        rational_field = FractionField(polynomial_ring)
        coordinate_ring = self.curve.coordinate_ring
        x_coordinate, y_coordinate = self.curve.coordinate_gens

        result = 0 * self.curve.dx
        multiplicative_order = Integers(self.curve.exponent)(characteristic).multiplicative_order()
        multiplier = (characteristic ** multiplicative_order - 1) // self.curve.exponent
        curve_polynomial = polynomial_ring(self.curve.polynomial)

        for denominator_power in range(0, int(self.curve.exponent)):
            component = _coerce_univariate_rational(self.jth_component(denominator_power), rational_field)
            h = rational_field(component * curve_polynomial ** (multiplier * denominator_power))
            denominator = h.denominator()
            h = polynomial_ring(h * denominator ** characteristic)

            next_denominator_power = (
                characteristic ** (multiplicative_order - 1) * denominator_power
            ) % self.curve.exponent
            polynomial_power = floor(
                characteristic ** (multiplicative_order - 1) * denominator_power / self.curve.exponent
            )
            cartier_polynomial = polynomial_part(characteristic, h)
            if hasattr(base_field, "cardinality") and base_field.cardinality() != characteristic:
                cartier_polynomial = polynomial_ring(
                    sum(
                        cartier_polynomial[degree].nth_root(characteristic) * x_variable ** degree
                        for degree in range(0, max(cartier_polynomial.degree(), -1) + 1)
                    )
                )

            denominator_xy = (
                coordinate_ring(curve_polynomial(x_coordinate)) ** polynomial_power
                * y_coordinate ** next_denominator_power
                * coordinate_ring(denominator(x_coordinate))
            )
            result += SuperellipticForm(
                self.curve,
                self.curve.function_field(cartier_polynomial(x_coordinate))
                / self.curve.function_field(denominator_xy),
            )
        return result

    def serre_duality_pairing(self, function, prec=20):
        '''Pair this form with a class in H^1(X, O_X) by residues at infinity.'''
        result = self.curve.base_ring.zero()
        for place in range(int(self.curve.nb_of_pts_at_infty)):
            result += (function * self).expansion_at_infty(place=place, prec=prec)[-1]
        return -result

    def coordinates(self, basis=0):
        '''Return coordinates of this holomorphic form in a given basis.'''
        if basis == 0:
            basis = self.curve.holomorphic_differentials_basis()
        _function_field, coordinate_ring, _x_coordinate, _y_coordinate = self.curve.fct_field
        common_denominator = LCM([omega.form.denominator() for omega in basis])
        basis_without_denominators = [common_denominator * omega.form for omega in basis]
        form_without_denominator = common_denominator * self.form
        return linear_representation_polynomials(
            coordinate_ring(form_without_denominator),
            [coordinate_ring(omega) for omega in basis_without_denominators],
        )

    def jth_component(self, denominator_power):
        '''Return h_j(x) when this form is written as sum_j h_j(x) dx / y^j.'''
        m = int(self.curve.exponent)
        _function_field, _coordinate_ring, x_coordinate, y_coordinate = self.curve.fct_field
        reduced = reduction(self.curve, y_coordinate ** m * self.form)
        polynomial_y, _polynomial_ring_y, x_auxiliary, _y_auxiliary, _x, _y = (
            _expression_to_polynomial_in_y(self.curve, reduced)
        )
        rational_field = polynomial_y.parent().base_ring()
        if denominator_power == 0:
            return polynomial_y[0] / rational_field(self.curve.polynomial(x_auxiliary))
        return polynomial_y[m - int(denominator_power)]

    def is_regular_on_U0(self):
        '''Return whether this form is regular on the affine chart U0.'''
        for denominator_power in range(0, int(self.curve.exponent)):
            component = self.jth_component(denominator_power)
            if component.denominator().degree() > 0:
                return False
        return True

    def is_regular_on_Uinfty(self):
        '''Return whether this form is regular on the chart containing infinity.'''
        parameters = self.curve.infinity_parameters()
        for denominator_power in range(1, int(self.curve.exponent)):
            component = self.jth_component(denominator_power)
            degree = degree_of_rational_function(component, self.curve.base_ring)
            if -degree * parameters.x_pole_order + denominator_power * parameters.y_pole_order - (
                parameters.x_pole_order + 1
            ) < 0:
                return False
        return True

    def expansion_at_infty(self, place=0, prec=10):
        '''Return the Laurent expansion at a place above infinity.'''
        from ascovers.superelliptic.function import SuperellipticFunction

        if not self.curve.dx_series or prec > self.curve.precision:
            self.curve.compute_expansions_at_infinity(prec)
        series_ring = LaurentSeriesRing(self.curve.base_ring, "t", default_prec=int(prec))
        coefficient = SuperellipticFunction(self.curve, self.form)
        coefficient_series = series_ring(coefficient.expansion_at_infty(place=place, prec=prec))
        return coefficient_series * series_ring(self.curve.dx_series[int(place)])

    def expansion(self, point, prec=50):
        '''Return the Laurent expansion in the completed local ring at a point.'''
        from ascovers.superelliptic.function import SuperellipticFunction

        if point in ZZ:
            return self.expansion_at_infty(place=point, prec=prec)
        dx_series = self.curve.x.expansion(point=point, prec=prec).derivative()
        coefficient = SuperellipticFunction(self.curve, self.form)
        return coefficient.expansion(point=point, prec=prec) * dx_series

    def residue(self, place=0, prec=30):
        '''Return the residue at a place above infinity.'''
        return self.expansion_at_infty(place=place, prec=prec)[-1]

    def reduce(self):
        '''Return this form with its coefficient reduced again.'''
        return SuperellipticForm(self.curve, reduction(self.curve, self.form))

    def reduce2(self):
        '''Return this form after reducing y^m times its coefficient.'''
        _function_field, _coordinate_ring, _x_coordinate, y_coordinate = self.curve.fct_field
        expression = reduction(self.curve, y_coordinate ** self.curve.exponent * self.form)
        return SuperellipticForm(self.curve, expression / y_coordinate ** self.curve.exponent)

    def integral(self):
        '''Compute a formal integral of this form in characteristic p.'''
        from ascovers.superelliptic.function import SuperellipticFunction

        characteristic = self.curve.characteristic
        if characteristic == 0:
            raise ValueError("this integral formula is only implemented in positive characteristic")

        base_field = self.curve.base_ring
        polynomial_ring = PolynomialRing(base_field, "x")
        rational_field = FractionField(polynomial_ring)
        coordinate_ring = self.curve.coordinate_ring
        x_coordinate, y_coordinate = self.curve.coordinate_gens
        curve_polynomial = polynomial_ring(self.curve.polynomial)

        result = 0 * self.curve.x
        multiplicative_order = Integers(self.curve.exponent)(characteristic).multiplicative_order()
        multiplier = (characteristic ** multiplicative_order - 1) // self.curve.exponent

        for denominator_power in range(0, int(self.curve.exponent)):
            component = _coerce_univariate_rational(self.jth_component(denominator_power), rational_field)
            h = rational_field(component * curve_polynomial ** (multiplier * denominator_power))
            denominator = h.denominator()
            h = polynomial_ring(h * denominator ** characteristic)

            next_denominator_power = (
                characteristic ** multiplicative_order * denominator_power
            ) % self.curve.exponent
            polynomial_power = floor(
                characteristic ** multiplicative_order * denominator_power / self.curve.exponent
            )
            denominator_xy = (
                coordinate_ring(curve_polynomial(x_coordinate)) ** polynomial_power
                * y_coordinate ** next_denominator_power
                * coordinate_ring(denominator(x_coordinate)) ** characteristic
            )
            result += SuperellipticFunction(
                self.curve,
                self.curve.function_field(h.integral()(x_coordinate))
                / self.curve.function_field(denominator_xy),
            )
        return result

    def inv_cartier(self):
        '''Return a preimage under Cartier for regular forms once regular forms are migrated.'''
        if not hasattr(self, "regular_form"):
            raise PendingMigrationError("regular forms have not been migrated yet")
        regular_form = self.regular_form()
        characteristic = self.curve.characteristic
        return (
            regular_form.dx ** characteristic * self.curve.x ** (characteristic - 1) * self.curve.dx
            + regular_form.dy ** characteristic * self.curve.y ** (characteristic - 1) * self.curve.y.diffn()
        )

    def valuation(self, place=0):
        '''Return the valuation at a place above infinity.'''
        series_ring = LaurentSeriesRing(self.curve.base_ring, "t")
        return series_ring(self.expansion_at_infty(place=place)).valuation()


superelliptic_form = SuperellipticForm
