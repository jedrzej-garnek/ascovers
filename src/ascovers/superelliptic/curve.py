'''Superelliptic curves and reductions modulo the equation y^m = f(x).'''

from dataclasses import dataclass
from importlib import import_module

from sage.all import (
    FractionField,
    HyperellipticCurve,
    LaurentSeriesRing,
    PolynomialRing,
    ZZ,
    gcd,
    matrix,
    vector,
    xgcd,
)

from ascovers.linear_algebra import flag
from ascovers.superelliptic.exceptions import PendingMigrationError
from ascovers.superelliptic.utils import cut, naive_hensel, root_of_unity, roots_in_base_field


@dataclass(frozen=True)
class InfinityParameters:
    '''Combinatorial data for the coordinate chart at infinity.'''

    number_of_places: object
    x_uniformizer_exponent: object
    y_uniformizer_exponent: object
    x_pole_order: object
    y_pole_order: object


class SuperellipticCurve:
    '''A superelliptic curve with affine equation y^m = f(x).'''

    def __init__(self, polynomial, exponent, prec=100, compute_expansions=True, variable_name="x"):
        '''Create the superelliptic curve y^m = f(x) over the coefficient field of f.

        INPUT:

        - ``polynomial`` -- a univariate Sage polynomial over a field.
        - ``exponent`` -- the positive integer ``m`` in the equation ``y^m = f(x)``.
        - ``prec`` -- precision used for Laurent expansions at infinity.
        - ``compute_expansions`` -- whether to compute infinity expansions immediately.
        - ``variable_name`` -- name of the normalized polynomial variable.

        EXAMPLES::

            sage: from ascovers import SuperellipticCurve
            sage: F = GF(3)
            sage: R.<x> = PolynomialRing(F)
            sage: C = SuperellipticCurve(x^3 - x, 2)
            sage: C.genus()
            1
        '''
        exponent = ZZ(exponent)
        if exponent <= 0:
            raise ValueError("the exponent must be a positive integer")

        input_parent = polynomial.parent()
        if not hasattr(input_parent, "base_ring"):
            raise TypeError("the defining equation must be a univariate polynomial over a field")
        coefficient_field = input_parent.base_ring()
        if hasattr(coefficient_field, "is_field") and not coefficient_field.is_field():
            raise TypeError("the defining polynomial must have coefficients in a field")

        characteristic = coefficient_field.characteristic()
        if characteristic != 0 and gcd(characteristic, exponent) != 1:
            raise ValueError("the exponent must be prime to the characteristic of the base field")

        self.polynomial_ring = PolynomialRing(coefficient_field, variable_name)
        self.variable = self.polynomial_ring.gen()
        self.polynomial = self.polynomial_ring(polynomial)
        if self.polynomial.degree() <= 0:
            raise ValueError("the defining polynomial must have positive degree")

        self.exponent = exponent
        self.degree = self.polynomial.degree()
        self._base_ring = coefficient_field
        self.precision = int(prec)
        self.prec = self.precision
        self.characteristic = characteristic
        self.nb_of_pts_at_infty = gcd(self.degree, self.exponent)

        self.coordinate_ring = PolynomialRing(coefficient_field, (variable_name, "y"))
        self.coordinate_x, self.coordinate_y = self.coordinate_ring.gens()
        self.function_field = FractionField(self.coordinate_ring)
        self.fct_field = (self.function_field, self.coordinate_ring, self.coordinate_x, self.coordinate_y)

        self.x_series = []
        self.y_series = []
        self.dx_series = []
        if compute_expansions:
            self.compute_expansions_at_infinity(self.precision)

    @property
    def base_ring(self):
        '''Return the coefficient field of the defining polynomial.'''
        return self._base_ring

    @property
    def x(self):
        '''Return the coordinate function x on the curve.'''
        return self.function(self.coordinate_x)

    @property
    def y(self):
        '''Return the coordinate function y on the curve.'''
        return self.function(self.coordinate_y)

    @property
    def dx(self):
        '''Return the differential dx on the curve.'''
        return self.form(self.coordinate_ring.one())

    @property
    def one(self):
        '''Return the constant function 1 on the curve.'''
        return self.function(self.coordinate_ring.one())

    @property
    def coordinate_gens(self):
        '''Return the coordinate-ring generators x and y.'''
        return self.coordinate_x, self.coordinate_y

    def __repr__(self):
        '''Return the standard textual representation of the curve.'''
        return (
            "Superelliptic curve with the equation y^"
            + str(self.exponent)
            + " = "
            + str(self.polynomial)
            + " over "
            + str(self.base_ring)
        )

    def _load_migrated_class(self, module_names, class_names, description):
        '''Load a migrated class that is expected in a later package module.'''
        for module_name in module_names:
            qualified_name = f"{__package__}.{module_name}"
            try:
                module = import_module(qualified_name)
            except ModuleNotFoundError as error:
                if error.name == qualified_name:
                    continue
                raise
            for class_name in class_names:
                if hasattr(module, class_name):
                    return getattr(module, class_name)

        expected = ", ".join(f"{module}.{class_name}" for module in module_names for class_name in class_names)
        raise PendingMigrationError(f"{description} could not be loaded; expected one of: {expected}")

    def _function_class(self):
        '''Return the migrated function class for this curve.'''
        return self._load_migrated_class(
            ("function", "functions", "superelliptic_function"),
            ("SuperellipticFunction", "superelliptic_function"),
            "Functions on superelliptic curves",
        )

    def _form_class(self):
        '''Return the migrated differential-form class for this curve.'''
        return self._load_migrated_class(
            ("form", "forms", "superelliptic_form"),
            ("SuperellipticForm", "superelliptic_form"),
            "Differential forms on superelliptic curves",
        )

    def _cech_class(self):
        '''Return the migrated Cech-cocycle class for this curve.'''
        return self._load_migrated_class(
            ("cech", "superelliptic_cech"),
            (
                "SuperellipticDeRhamCocycle",
                "SuperellipticCechCocycle",
                "SuperellipticCech",
                "superelliptic_cech",
            ),
            "Cech cocycles on superelliptic curves",
        )

    def function(self, expression):
        '''Coerce an expression to a function on this superelliptic curve.'''
        function_class = self._function_class()
        return function_class(self, expression)

    def form(self, expression):
        '''Coerce an expression to a differential form on this superelliptic curve.'''
        form_class = self._form_class()
        return form_class(self, expression)

    def cech(self, omega, function):
        '''Build a Cech representative for de Rham cohomology.'''
        cech_class = self._cech_class()
        return cech_class(self, omega, function)

    def infinity_parameters(self):
        '''Return normalized Bezout data used for expansions at infinity.'''
        number_of_places, bezout_exponent, bezout_degree = xgcd(self.exponent, self.degree)
        x_uniformizer_exponent = -bezout_exponent
        y_uniformizer_exponent = bezout_degree
        x_pole_order = self.exponent // number_of_places
        y_pole_order = self.degree // number_of_places

        while x_uniformizer_exponent < 0:
            x_uniformizer_exponent += y_pole_order
            y_uniformizer_exponent += x_pole_order

        return InfinityParameters(
            number_of_places=number_of_places,
            x_uniformizer_exponent=x_uniformizer_exponent,
            y_uniformizer_exponent=y_uniformizer_exponent,
            x_pole_order=x_pole_order,
            y_pole_order=y_pole_order,
        )

    def reciprocal_polynomial(self):
        '''Return x^r f(1/x), where r is the degree of f.'''
        fraction_field = FractionField(self.polynomial_ring)
        variable = fraction_field(self.variable)
        reciprocal = variable ** self.degree * self.polynomial(1 / variable)
        return self.polynomial_ring(reciprocal)

    def _infinity_start_values(self, parameters):
        '''Return initial roots for Hensel lifts at the places above infinity.'''
        leading_coefficient = self.polynomial.leading_coefficient()
        if leading_coefficient == self.base_ring.one():
            primitive_root = root_of_unity(self.base_ring, parameters.number_of_places)
            return [primitive_root ** place for place in range(int(parameters.number_of_places))]
        return roots_in_base_field(self.base_ring, parameters.number_of_places, leading_coefficient)

    def compute_expansions_at_infinity(self, precision=None):
        '''Compute Laurent expansions of x, y, and dx at every rational place at infinity.'''
        if precision is None:
            precision = self.precision
        precision = int(precision)

        parameters = self.infinity_parameters()
        series_ring = LaurentSeriesRing(self.base_ring, "t", default_prec=precision)
        series_variable = series_ring.gen()
        polynomial_ring = PolynomialRing(series_ring, "W")
        polynomial_variable = polynomial_ring.gen()
        polynomial_field = FractionField(polynomial_ring)

        reciprocal = self.reciprocal_polynomial()
        transformed_polynomial = (
            polynomial_field(
                reciprocal(series_variable ** parameters.x_pole_order
                           * polynomial_variable ** parameters.y_uniformizer_exponent)
            )
            - polynomial_variable ** parameters.number_of_places
        )

        self.x_series = []
        self.y_series = []
        self.dx_series = []
        for start_value in self._infinity_start_values(parameters):
            lifted_root = naive_hensel(transformed_polynomial, self.base_ring, start=start_value, prec=precision)
            x_series = series_ring(
                1 / (series_variable ** parameters.x_pole_order
                     * lifted_root ** parameters.y_uniformizer_exponent)
            )
            y_series = series_ring(
                1 / (series_variable ** parameters.y_pole_order
                     * lifted_root ** parameters.x_uniformizer_exponent)
            )
            self.x_series.append(x_series)
            self.y_series.append(y_series)
            self.dx_series.append(x_series.derivative())

        return self.x_series, self.y_series, self.dx_series

    def _holomorphic_differential_degrees(self):
        '''Return pairs (i, j) for basis elements x^i dx / y^j.'''
        degrees = []
        number_of_places = self.nb_of_pts_at_infty
        for denominator_power in range(1, int(self.exponent)):
            for shifted_x_degree in range(1, int(self.degree)):
                if self.degree * denominator_power - self.exponent * shifted_x_degree >= number_of_places:
                    degrees.append((shifted_x_degree - 1, denominator_power))
        return degrees

    def basis_holomorphic_differentials_degree(self):
        '''Return the holomorphic differential basis and its degree dictionary.'''
        function_field, _coordinate_ring, x_coordinate, y_coordinate = self.fct_field
        basis = [
            self.form(function_field(x_coordinate ** x_degree / y_coordinate ** denominator_power))
            for x_degree, denominator_power in self._holomorphic_differential_degrees()
        ]
        degrees = {
            index: degree_pair
            for index, degree_pair in enumerate(self._holomorphic_differential_degrees())
        }
        return basis, degrees

    def holomorphic_differentials_basis(self):
        '''Return a basis of holomorphic differentials on the curve.'''
        basis, _degrees = self.basis_holomorphic_differentials_degree()
        return basis

    def degrees_holomorphic_differentials(self):
        '''Return degree pairs (i, j) such that x^i dx / y^j is holomorphic.'''
        return {
            index: degree_pair
            for index, degree_pair in enumerate(self._holomorphic_differential_degrees())
        }

    def basis_de_rham_degrees(self):
        '''Return a de Rham basis together with degree data on the two standard opens.'''
        genus = int(self.genus())
        basis_holomorphic = self.holomorphic_differentials_basis()
        basis = [None] * (2 * genus)

        zero_function = self.function(self.function_field.zero())
        for index, differential in enumerate(basis_holomorphic):
            basis[index] = self.cech(differential, zero_function)

        degrees_holomorphic = self.degrees_holomorphic_differentials()
        inverse_degrees = {degree: index for index, degree in degrees_holomorphic.items()}
        degrees_on_affine = {}
        degrees_at_infinity = {}

        function_field, _coordinate_ring, x_coordinate, y_coordinate = self.fct_field
        polynomial_variable = self.variable

        for denominator_power in range(1, int(self.exponent)):
            for shifted_x_degree in range(1, int(self.degree)):
                condition = (
                    self.degree * (self.exponent - denominator_power)
                    - self.exponent * shifted_x_degree
                    >= self.nb_of_pts_at_infty
                )
                if not condition:
                    continue

                holomorphic_index = inverse_degrees[(shifted_x_degree - 1, self.exponent - denominator_power)]
                cohomology_function = self.function(
                    function_field(
                        self.exponent
                        * y_coordinate ** (self.exponent - denominator_power)
                        / x_coordinate ** shifted_x_degree
                    )
                )
                pairing_constant = basis_holomorphic[holomorphic_index].serre_duality_pairing(cohomology_function)

                relation_polynomial = (
                    (self.exponent - denominator_power) * polynomial_variable * self.polynomial.derivative()
                    - self.exponent * shifted_x_degree * self.polynomial
                )
                psi = cut(relation_polynomial, shifted_x_degree)
                affine_form = self.form(function_field(psi(x_coordinate) / y_coordinate ** denominator_power))
                basis_index = holomorphic_index + genus
                basis[basis_index] = (1 / pairing_constant) * self.cech(affine_form, cohomology_function)
                degrees_on_affine[basis_index] = (psi.degree(), denominator_power)
                degrees_at_infinity[basis_index] = (-shifted_x_degree, self.exponent - denominator_power)

        if any(element is None for element in basis):
            raise ValueError("failed to construct a complete de Rham basis")
        return basis, degrees_on_affine, degrees_at_infinity

    def de_rham_basis(self):
        '''Return a basis of algebraic de Rham cohomology represented by Cech cocycles.'''
        basis, _degrees_on_affine, _degrees_at_infinity = self.basis_de_rham_degrees()
        return basis

    def degrees_de_rham0(self):
        '''Return degree data for the affine component of the de Rham basis.'''
        _basis, degrees_on_affine, _degrees_at_infinity = self.basis_de_rham_degrees()
        return degrees_on_affine

    def degrees_de_rham1(self):
        '''Return degree data for the infinity component of the de Rham basis.'''
        _basis, _degrees_on_affine, degrees_at_infinity = self.basis_de_rham_degrees()
        return degrees_at_infinity

    def is_smooth(self):
        '''Return whether the affine model has squarefree defining polynomial.'''
        return self.polynomial.discriminant() != 0

    def genus(self):
        '''Return the genus of the smooth projective superelliptic curve.'''
        genus_numerator = (
            (self.degree - 1) * (self.exponent - 1)
            - self.nb_of_pts_at_infty
            + 1
        )
        return ZZ(genus_numerator // 2)

    def verschiebung_matrix(self):
        '''Return the Verschiebung matrix on de Rham cohomology in the standard basis.'''
        basis = self.de_rham_basis()
        genus = int(self.genus())
        matrix_result = matrix(self.base_ring, 2 * genus, 2 * genus)
        for row_index, cohomology_class in enumerate(basis):
            matrix_result[row_index, :] = cohomology_class.verschiebung().coordinates()
        return matrix_result

    def dr_frobenius_matrix(self):
        '''Return the Frobenius matrix on de Rham cohomology in the standard basis.'''
        basis = self.de_rham_basis()
        genus = int(self.genus())
        matrix_result = matrix(self.base_ring, 2 * genus, 2 * genus)
        for row_index, cohomology_class in enumerate(basis):
            matrix_result[row_index, :] = cohomology_class.frobenius().coordinates()
        return matrix_result

    def cartier_matrix(self):
        '''Return the Cartier matrix on holomorphic differentials.'''
        basis = self.holomorphic_differentials_basis()
        genus = int(self.genus())
        matrix_result = matrix(self.base_ring, genus, genus)
        for column_index, differential in enumerate(basis):
            matrix_result[:, column_index] = vector(differential.cartier().coordinates())
        return matrix_result

    def frobenius_matrix(self, prec=20):
        '''Return the Frobenius matrix on H^1(X, O_X).'''
        genus = int(self.genus())
        characteristic = self.base_ring.characteristic()
        matrix_result = matrix(self.base_ring, genus, genus)
        for row_index, function in enumerate(self.cohomology_of_structure_sheaf_basis()):
            matrix_result[row_index, :] = vector((function ** characteristic).coordinates(prec=prec))
        return matrix_result.transpose()

    def p_rank(self):
        '''Return the p-rank for hyperelliptic curves, where Sage has a native routine.'''
        if self.exponent != 2:
            raise NotImplementedError("p-rank is currently implemented only for hyperelliptic curves")
        polynomial_ring = PolynomialRing(self.base_ring, "t")
        hyperelliptic_polynomial = polynomial_ring(self.polynomial)
        hyperelliptic_curve = HyperellipticCurve(hyperelliptic_polynomial, 0)
        return hyperelliptic_curve.p_rank()

    def a_number(self):
        '''Return the a-number computed from the rank of the Cartier matrix.'''
        return self.genus() - self.cartier_matrix().rank()

    def final_type(self, test=False):
        '''Return the final type computed from Frobenius and Verschiebung.'''
        return flag(self.frobenius_matrix(), self.verschiebung_matrix(), self.characteristic, test)

    def cohomology_of_structure_sheaf_basis(self):
        '''Return the Serre-dual basis of H^1(X, O_X).'''
        genus = int(self.genus())
        basis = [None] * genus
        degrees = self.degrees_holomorphic_differentials()
        inverse_degrees = {degree: index for index, degree in degrees.items()}
        basis_holomorphic = self.holomorphic_differentials_basis()
        function_field, _coordinate_ring, x_coordinate, y_coordinate = self.fct_field

        for denominator_power in range(1, int(self.exponent)):
            for shifted_x_degree in range(1, int(self.degree)):
                condition = (
                    self.degree * (self.exponent - denominator_power)
                    - self.exponent * shifted_x_degree
                    >= self.nb_of_pts_at_infty
                )
                if not condition:
                    continue

                index = inverse_degrees[(shifted_x_degree - 1, self.exponent - denominator_power)]
                cohomology_function = self.function(
                    function_field(
                        self.exponent
                        * y_coordinate ** (self.exponent - denominator_power)
                        / x_coordinate ** shifted_x_degree
                    )
                )
                pairing_constant = basis_holomorphic[index].serre_duality_pairing(cohomology_function)
                basis[index] = (1 / pairing_constant) * cohomology_function

        if any(element is None for element in basis):
            raise ValueError("failed to construct a complete H^1(X, O_X) basis")
        return basis

    def uniformizer(self, place=0):
        '''Return a function with valuation one at a place above infinity.'''
        parameters = self.infinity_parameters()
        uniformizer = self.one
        if parameters.x_uniformizer_exponent:
            uniformizer *= self.x ** int(parameters.x_uniformizer_exponent)
        if parameters.y_uniformizer_exponent:
            uniformizer /= self.y ** int(parameters.y_uniformizer_exponent)
        return uniformizer

    def _load_holomorphic_combination_helpers(self):
        '''Load Riemann-Roch helper functions shared with AS covers.'''
        from ascovers.as_covers.holomorphic import (
            holomorphic_combinations_fcts,
            holomorphic_combinations_forms,
        )

        return holomorphic_combinations_fcts, holomorphic_combinations_forms

    def riemann_roch_space(self, pole_orders, threshold=8):
        '''Find functions with bounded pole orders at the places above infinity.'''
        holomorphic_combinations_fcts, _holomorphic_combinations_forms = self._load_holomorphic_combination_helpers()
        function_field, _coordinate_ring, x_coordinate, y_coordinate = self.fct_field

        candidates = []
        for x_power in range(0, threshold * int(self.degree)):
            for y_power in range(0, int(self.exponent)):
                candidate = self.function(function_field(x_coordinate ** x_power * y_coordinate ** y_power))
                candidates.append((candidate, candidate.expansion_at_infty()))

        functions = holomorphic_combinations_fcts(candidates, pole_orders[0])
        for place in range(1, int(self.nb_of_pts_at_infty)):
            functions = [(function, function.expansion_at_infty(place=place)) for function in functions]
            functions = holomorphic_combinations_fcts(functions, pole_orders[place])
        return functions

    def riemann_roch_space_forms(self, pole_orders, threshold=8):
        '''Find differential forms with bounded pole orders at the places above infinity.'''
        _holomorphic_combinations_fcts, holomorphic_combinations_forms = self._load_holomorphic_combination_helpers()
        function_field, _coordinate_ring, x_coordinate, y_coordinate = self.fct_field

        candidates = []
        for x_power in range(0, threshold * int(self.degree)):
            for denominator_power in range(0, int(self.exponent)):
                candidate = self.form(function_field(x_coordinate ** x_power / y_coordinate ** denominator_power))
                candidates.append((candidate, candidate.expansion_at_infty()))

        forms = holomorphic_combinations_forms(candidates, pole_orders[0])
        for place in range(1, int(self.nb_of_pts_at_infty)):
            forms = [(form, form.expansion_at_infty(place=place)) for form in forms]
            forms = holomorphic_combinations_forms(forms, pole_orders[place])
        return forms


def _rational_x_to_function_field(coefficient, function_field, x_coordinate):
    '''Evaluate a rational function in one variable at the coordinate x.'''
    numerator = coefficient.numerator()
    denominator = coefficient.denominator()
    return function_field(numerator(x_coordinate)) / function_field(denominator(x_coordinate))


def _polynomial_y_to_function_field(polynomial, function_field, x_coordinate, y_coordinate):
    '''Evaluate a polynomial in y with rational x-coefficients in the curve function field.'''
    result = function_field.zero()
    for degree, coefficient in polynomial.dict().items():
        y_degree = degree[0] if isinstance(degree, tuple) else degree
        result += (
            _rational_x_to_function_field(coefficient, function_field, x_coordinate)
            * function_field(y_coordinate ** int(y_degree))
        )
    return function_field(result)


def _expression_to_polynomial_in_y(curve, expression):
    '''Convert a reduced function-field expression to a polynomial in y over F(x).'''
    function_field, coordinate_ring, x_coordinate, y_coordinate = curve.fct_field
    fraction = function_field(expression)
    numerator = coordinate_ring(fraction.numerator())
    denominator = coordinate_ring(fraction.denominator())

    polynomial_ring_x = PolynomialRing(curve.base_ring, "x1")
    x_auxiliary = polynomial_ring_x.gen()
    rational_field_x = FractionField(polynomial_ring_x)
    polynomial_ring_y = PolynomialRing(rational_field_x, "y1")
    y_auxiliary = polynomial_ring_y.gen()
    polynomial_field_y = FractionField(polynomial_ring_y)

    numerator_y = polynomial_ring_y(numerator(x_auxiliary, y_auxiliary))
    denominator_y = polynomial_ring_y(denominator(x_auxiliary, y_auxiliary))
    polynomial_y = polynomial_ring_y(polynomial_field_y(numerator_y) / polynomial_field_y(denominator_y))
    return polynomial_y, polynomial_ring_y, x_auxiliary, y_auxiliary, x_coordinate, y_coordinate


def reduction(curve, expression):
    '''Reduce a rational expression modulo y^m = f(x) to degree less than m in y.'''
    function_field, coordinate_ring, x_coordinate, y_coordinate = curve.fct_field
    fraction = function_field(expression)
    numerator = coordinate_ring(fraction.numerator())
    denominator = coordinate_ring(fraction.denominator())

    polynomial_ring_x = PolynomialRing(curve.base_ring, "x1")
    x_auxiliary = polynomial_ring_x.gen()
    rational_field_x = FractionField(polynomial_ring_x)
    polynomial_ring_y = PolynomialRing(rational_field_x, "y1")
    y_auxiliary = polynomial_ring_y.gen()
    polynomial_field_y = FractionField(polynomial_ring_y)

    numerator_y = polynomial_ring_y(numerator(x_auxiliary, y_auxiliary))
    denominator_y = polynomial_ring_y(denominator(x_auxiliary, y_auxiliary))
    relation = y_auxiliary ** curve.exponent - polynomial_ring_y(curve.polynomial(x_auxiliary))
    gcd_denominator, inverse_denominator, _relation_coefficient = xgcd(denominator_y, relation)
    reduced_polynomial = polynomial_ring_y(
        polynomial_field_y(numerator_y) * inverse_denominator / gcd_denominator
    )

    while reduced_polynomial.degree() >= curve.exponent:
        current_degree = reduced_polynomial.degree()
        leading_component = reduced_polynomial[current_degree]
        quotient = current_degree // curve.exponent
        remainder = current_degree % curve.exponent
        reduced_polynomial = polynomial_ring_y(
            reduced_polynomial
            - leading_component * y_auxiliary ** current_degree
            + curve.polynomial(x_auxiliary) ** quotient * leading_component * y_auxiliary ** remainder
        )

    return _polynomial_y_to_function_field(reduced_polynomial, function_field, x_coordinate, y_coordinate)


def reduction_form(curve, expression):
    '''Reduce a differential-form coefficient to the standard basis with denominators y^j.'''
    function_field, _coordinate_ring, x_coordinate, y_coordinate = curve.fct_field
    reduced_expression = reduction(curve, expression)
    reduced_polynomial, _polynomial_ring_y, _x_auxiliary, _y_auxiliary, _x, _y = (
        _expression_to_polynomial_in_y(curve, reduced_expression)
    )

    curve_polynomial = function_field(curve.polynomial(x_coordinate))
    result = function_field.zero()
    for y_power in range(0, int(curve.exponent)):
        coefficient = reduced_polynomial[y_power]
        if coefficient == 0:
            continue
        coefficient_xy = _rational_x_to_function_field(coefficient, function_field, x_coordinate)
        if y_power == 0:
            result += coefficient_xy
        else:
            result += coefficient_xy * curve_polynomial / function_field(y_coordinate ** (curve.exponent - y_power))
    return function_field(result)


superelliptic = SuperellipticCurve
