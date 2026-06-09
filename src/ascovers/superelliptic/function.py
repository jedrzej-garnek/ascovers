'''Rational functions on superelliptic curves.'''

from sage.all import LaurentSeriesRing, PolynomialRing, ZZ

from ascovers.superelliptic.curve import _expression_to_polynomial_in_y, reduction
from ascovers.superelliptic.utils import naive_hensel, reverse_power_series


class SuperellipticFunction:
    '''A rational function on a superelliptic curve.'''

    def __init__(self, curve, expression):
        '''Create the function represented by an expression in x and y.'''
        if isinstance(expression, SuperellipticFunction):
            if expression.curve is not curve:
                raise TypeError("cannot coerce a function from a different curve")
            expression = expression.function
        self.curve = curve
        self.function = reduction(curve, expression)

    def __eq__(self, other):
        '''Return whether two superelliptic functions have the same reduced expression.'''
        if not isinstance(other, SuperellipticFunction):
            try:
                other = SuperellipticFunction(self.curve, other)
            except (TypeError, ValueError):
                return False
        return self.curve is other.curve and self.function == other.function

    def __repr__(self):
        '''Return the textual representation of the reduced expression.'''
        return str(self.function)

    def _coerce_expression(self, other):
        '''Return the raw expression underlying another function or scalar.'''
        if isinstance(other, SuperellipticFunction):
            if other.curve is not self.curve:
                raise TypeError("cannot combine functions on different curves")
            return other.function
        return other

    def jth_component(self, exponent):
        '''Return the coefficient of y^exponent in the reduced function.'''
        polynomial_y, _polynomial_ring_y, _x_auxiliary, _y_auxiliary, _x, _y = (
            _expression_to_polynomial_in_y(self.curve, self.function)
        )
        return polynomial_y[int(exponent)]

    def __add__(self, other):
        '''Add two functions or a function and a scalar expression.'''
        return SuperellipticFunction(self.curve, self.function + self._coerce_expression(other))

    def __radd__(self, other):
        '''Add a scalar expression and this function.'''
        return self + other

    def __neg__(self):
        '''Return the additive inverse of this function.'''
        return SuperellipticFunction(self.curve, -self.function)

    def __sub__(self, other):
        '''Subtract another function or scalar expression.'''
        return SuperellipticFunction(self.curve, self.function - self._coerce_expression(other))

    def __rsub__(self, other):
        '''Subtract this function from a scalar expression.'''
        return SuperellipticFunction(self.curve, self._coerce_expression(other) - self.function)

    def __mul__(self, other):
        '''Multiply by a function, scalar expression, or differential form.'''
        from ascovers.superelliptic.form import SuperellipticForm

        if isinstance(other, SuperellipticForm):
            if other.curve is not self.curve:
                raise TypeError("cannot multiply objects on different curves")
            return SuperellipticForm(self.curve, self * other.form)
        return SuperellipticFunction(self.curve, self.function * self._coerce_expression(other))

    def __rmul__(self, other):
        '''Multiply this function by a scalar expression.'''
        return self * other

    def __truediv__(self, other):
        '''Divide by another function or scalar expression.'''
        return SuperellipticFunction(self.curve, self.function / self._coerce_expression(other))

    def __rtruediv__(self, other):
        '''Divide a scalar expression by this function.'''
        return SuperellipticFunction(self.curve, self._coerce_expression(other) / self.function)

    def __pow__(self, exponent):
        '''Raise this function to an integer power.'''
        return SuperellipticFunction(self.curve, self.function ** exponent)

    def diffn(self):
        '''Return the differential of this function as a multiple of dx.'''
        from ascovers.superelliptic.form import SuperellipticForm

        function_field, _coordinate_ring, x_coordinate, y_coordinate = self.curve.fct_field
        function = function_field(self.function)
        polynomial_derivative = self.curve.polynomial.derivative()(x_coordinate)
        derivative_x = function.derivative(x_coordinate)
        derivative_y = function.derivative(y_coordinate) * polynomial_derivative
        derivative_y /= self.curve.exponent * y_coordinate ** (self.curve.exponent - 1)
        return SuperellipticForm(self.curve, derivative_x + derivative_y)

    def coordinates(self, basis=0, prec=50):
        '''Return coordinates in H^1(X, O_X) dual to a basis of holomorphic forms.'''
        if basis == 0:
            basis = self.curve.holomorphic_differentials_basis()
        coordinates = [self.curve.base_ring.zero()] * int(self.curve.genus())
        for index, differential in enumerate(basis):
            coordinates[index] = differential.serre_duality_pairing(self, prec=prec)
        return coordinates

    def expansion_at_infty(self, place=0, prec=20):
        '''Return the Laurent expansion of this function at a place above infinity.'''
        place = int(place)
        if not self.curve.x_series or prec > self.curve.precision:
            self.curve.compute_expansions_at_infinity(prec)
        series_ring = LaurentSeriesRing(self.curve.base_ring, "t", default_prec=int(prec))
        x_series = series_ring(self.curve.x_series[place])
        y_series = series_ring(self.curve.y_series[place])
        return series_ring(self.function(x_series, y_series))

    def expansion(self, point, prec=50):
        '''Return the Laurent expansion in the completed local ring at a point.'''
        if point in ZZ:
            return self.expansion_at_infty(place=point, prec=prec)

        x0, y0 = point[0], point[1]
        series_ring = LaurentSeriesRing(self.curve.base_ring, "t", default_prec=int(prec))
        series_variable = series_ring.gen()
        polynomial_ring = PolynomialRing(series_ring, "W")
        polynomial_variable = polynomial_ring.gen()

        shifted_x = series_variable + x0
        derivative_at_x0 = self.curve.polynomial.derivative()(x0)
        if y0 != 0 or derivative_at_x0 == 0:
            y_series = naive_hensel(
                polynomial_variable ** self.curve.exponent - self.curve.polynomial(shifted_x),
                self.curve.base_ring,
                start=y0,
                prec=prec,
            )
            return series_ring(self.function(shifted_x, y_series))

        polynomial_variable_x = PolynomialRing(self.curve.base_ring, "x").gen()
        shifted_polynomial = self.curve.polynomial(polynomial_variable_x + x0)
        inverse_series = reverse_power_series(shifted_polynomial(series_variable), prec=prec)
        x_series = inverse_series(series_variable ** self.curve.exponent) + x0
        return series_ring(self.function(x_series, series_variable))

    def pth_root(self):
        '''Return the p-th root of this function when it is a p-th power.'''
        from ascovers.superelliptic.form import SuperellipticForm

        if self.diffn().form != 0:
            raise ValueError("function is not a p-th power")
        _function_field, _coordinate_ring, x_coordinate, _y_coordinate = self.curve.fct_field
        auxiliary_form = SuperellipticForm(self.curve, self.function / x_coordinate)
        auxiliary_form = auxiliary_form.cartier()
        auxiliary_form = self.curve.x * auxiliary_form
        return SuperellipticFunction(self.curve, auxiliary_form.form)

    def valuation(self, place=0):
        '''Return the valuation at a place above infinity.'''
        series_ring = LaurentSeriesRing(self.curve.base_ring, "t")
        return series_ring(self.expansion_at_infty(place=place)).valuation()

    def numerator(self):
        '''Return the numerator as a superelliptic function.'''
        function_field, _coordinate_ring, _x_coordinate, _y_coordinate = self.curve.fct_field
        return SuperellipticFunction(self.curve, function_field(self.function.numerator()))

    def denominator(self):
        '''Return the denominator as a superelliptic function.'''
        function_field, _coordinate_ring, _x_coordinate, _y_coordinate = self.curve.fct_field
        return SuperellipticFunction(self.curve, function_field(self.function.denominator()))

    def reduce(self):
        '''Return this function with its representative reduced again.'''
        return SuperellipticFunction(self.curve, self.function)

    def evaluate(self, point):
        '''Evaluate this function at an affine point when no cancellation is needed.'''
        if getattr(point, "is_at_infty", False):
            raise NotImplementedError("evaluation at infinity will be implemented with point objects")
        if hasattr(point, "xy"):
            x0, y0 = point.xy
        else:
            x0, y0 = point[0], point[1]
        numerator = self.function.numerator()
        denominator = self.function.denominator()
        denominator_value = denominator(x0, y0)
        if denominator_value == 0:
            raise ZeroDivisionError("the function has a pole or requires local cancellation at this point")
        return numerator(x0, y0) / denominator_value


superelliptic_function = SuperellipticFunction
