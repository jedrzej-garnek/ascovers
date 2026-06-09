'''Tests for superelliptic curves, functions, and forms.'''

import unittest

try:
    from sage.all import GF, HyperellipticCurve, PolynomialRing, vector
except ModuleNotFoundError:  # pragma: no cover - exercised outside Sage.
    GF = None
    PolynomialRing = None
    vector = None


if GF is not None:
    from ascovers import (
        SuperellipticCurve,
        SuperellipticForm,
        SuperellipticFunction,
        reduction,
        reduction_form,
        superelliptic,
    )


@unittest.skipIf(GF is None, "SageMath is required")
class SuperellipticCurveTest(unittest.TestCase):
    '''Tests for migrated superelliptic curve functionality.'''

    def test_basic_curve_invariants(self):
        '''Check the basic invariants of a smooth hyperelliptic example.'''
        field = GF(3)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 3 - x, 2, prec=20)

        self.assertEqual(curve.base_ring, field)
        self.assertEqual(curve.exponent, 2)
        self.assertEqual(curve.degree, 3)
        self.assertEqual(curve.nb_of_pts_at_infty, 1)
        self.assertEqual(curve.genus(), 1)
        self.assertTrue(curve.is_smooth())
        self.assertEqual(
            repr(curve),
            "Superelliptic curve with the equation y^2 = x^3 + 2*x over Finite Field of size 3",
        )

    def test_legacy_lowercase_alias(self):
        '''Check that the legacy constructor name remains available.'''
        field = GF(5)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = superelliptic(x ** 5 + x + 1, 4, compute_expansions=False)
        self.assertIsInstance(curve, SuperellipticCurve)

    def test_coordinate_function_objects_are_available(self):
        '''Check that curve coordinate shortcuts now create function and form objects.'''
        field = GF(3)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 3 - x, 2, compute_expansions=False)

        self.assertIsInstance(curve.x, SuperellipticFunction)
        self.assertIsInstance(curve.y, SuperellipticFunction)
        self.assertIsInstance(curve.dx, SuperellipticForm)

    def test_holomorphic_differential_basis_degrees(self):
        '''Check the holomorphic differential basis and degree dictionary.'''
        field = GF(7)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 5 + x, 4, compute_expansions=False)

        expected_degrees = {
            0: (0, 1),
            1: (0, 2),
            2: (1, 2),
            3: (0, 3),
            4: (1, 3),
            5: (2, 3),
        }
        self.assertEqual(curve.degrees_holomorphic_differentials(), expected_degrees)
        self.assertEqual(len(curve.holomorphic_differentials_basis()), curve.genus())

    def test_infinity_expansions_satisfy_curve_equation(self):
        '''Check that the computed Laurent expansions satisfy y^m = f(x).'''
        field = GF(3)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 3 - x, 2, prec=40)

        x_series = curve.x_series[0]
        y_series = curve.y_series[0]
        parameters = curve.infinity_parameters()

        self.assertEqual(y_series ** curve.exponent, curve.polynomial(x_series))
        self.assertEqual(x_series.valuation(), -parameters.x_pole_order)
        self.assertEqual(y_series.valuation(), -parameters.y_pole_order)

    def test_reduction_modulo_superelliptic_relation(self):
        '''Check reduction of functions modulo the equation y^m = f(x).'''
        field = GF(5)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 3 + x + 1, 2, compute_expansions=False)
        function_field, _coordinate_ring, x_coordinate, y_coordinate = curve.fct_field
        defining_polynomial = curve.polynomial(x_coordinate)

        reduced = reduction(curve, y_coordinate ** 3 + x_coordinate * y_coordinate ** 2)
        expected = function_field(defining_polynomial * y_coordinate + x_coordinate * defining_polynomial)

        self.assertEqual(reduced, expected)

    def test_reduction_form_standard_denominators(self):
        '''Check reduction of form coefficients to powers of y in the denominator.'''
        field = GF(5)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 3 + x + 1, 2, compute_expansions=False)
        function_field, _coordinate_ring, x_coordinate, y_coordinate = curve.fct_field
        defining_polynomial = curve.polynomial(x_coordinate)

        reduced = reduction_form(curve, y_coordinate)
        expected = function_field(defining_polynomial / y_coordinate)

        self.assertEqual(reduced, expected)

    def test_function_arithmetic_and_differential(self):
        '''Check arithmetic on functions and differentiation on the curve.'''
        field = GF(5)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 3 + x + 1, 2, compute_expansions=False)

        self.assertEqual(curve.y ** 2, curve.x ** 3 + curve.x + curve.one)
        self.assertEqual((curve.x + 2 * curve.y + curve.one).jth_component(1), 2)
        self.assertEqual(curve.x.diffn(), curve.dx)

    def test_form_coordinates(self):
        '''Check coordinates of a linear combination of the holomorphic basis.'''
        field = GF(7)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 5 + x, 4)
        basis = curve.holomorphic_differentials_basis()
        coefficients = [field(index + 1) for index in range(int(curve.genus()))]

        form = 0 * curve.dx
        for index in range(int(curve.genus())):
            form += coefficients[index] * basis[index]

        self.assertEqual(vector(form.coordinates(basis=basis)), vector(coefficients))

    def test_function_expansion_at_infinity(self):
        '''Check expansions of x and y at infinity.'''
        field = GF(3)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 3 - x, 2, prec=80)
        x_expansion = curve.x.expansion_at_infty(prec=80)
        y_expansion = curve.y.expansion_at_infty(prec=80)
        parameters = curve.infinity_parameters()

        self.assertEqual(y_expansion ** 2, curve.polynomial(x_expansion))
        self.assertEqual(x_expansion.valuation(), -parameters.x_pole_order)
        self.assertEqual(y_expansion.valuation(), -parameters.y_pole_order)

    def test_pth_root(self):
        '''Check p-th roots of functions when they exist.'''
        characteristic = 3
        field = GF(characteristic)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x ** 5 + x, 4)

        function = (curve.x ** 5) * (curve.y ** 2) + 2 * (curve.x ** 2) * (curve.y ** 3)
        self.assertEqual((function ** characteristic).pth_root(), function)
        with self.assertRaises(ValueError):
            curve.x.pth_root()

    def test_hyperelliptic_a_number_matches_sage(self):
        '''Check the a-number against Sage's native hyperelliptic implementation.'''
        field = GF(67)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        polynomial = x ** 7 + x ** 3 + x

        curve = SuperellipticCurve(polynomial, 2)
        sage_curve = HyperellipticCurve(polynomial)

        self.assertEqual(curve.a_number(), sage_curve.a_number())


if __name__ == "__main__":
    unittest.main()
