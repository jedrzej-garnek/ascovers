'''Tests for Artin-Schreier type covers.'''

import unittest

try:
    from sage.all import GF, PolynomialRing
except ModuleNotFoundError:  # pragma: no cover - exercised outside Sage.
    GF = None
    PolynomialRing = None


if GF is not None:
    from ascovers import (
        ArtinSchreierCover,
        ArtinSchreierDeRhamCocycle,
        ArtinSchreierForm,
        ArtinSchreierFunction,
        ArtinSchreierPolyform,
        CoverTemplate,
        SuperellipticCurve,
        as_cover,
        combination_components,
        dual_elt,
        elementary_cover,
        elementary_gp,
        elementary_template,
        heisenberg,
        magma_module_decomposition,
        heisenberg_template,
        quaternion_template,
        witt_template,
    )


@unittest.skipIf(GF is None, "SageMath is required")
class ArtinSchreierCoverTest(unittest.TestCase):
    '''Tests for migrated Artin-Schreier cover functionality.'''

    def quotient_curve(self):
        '''Return a small quotient curve used by the AS tests.'''
        field = GF(5)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        return SuperellipticCurve(x ** 3 + x + 1, 2, prec=60)

    def elementary_cover(self):
        '''Return a one-generator elementary cover for tests.'''
        curve = self.quotient_curve()
        return elementary_cover([curve.x], prec=60)

    def rational_elementary_cover(self):
        '''Return a small rational AS cover for basis-search tests.'''
        field = GF(2)
        polynomial_ring = PolynomialRing(field, "x")
        x = polynomial_ring.gen()
        curve = SuperellipticCurve(x, 1, prec=100)
        return elementary_cover([curve.x ** 3], prec=100)

    def test_group_arithmetic_and_legacy_aliases(self):
        '''Check finite group operations and old constructor aliases.'''
        group = elementary_gp(3, 2)
        first, second = group.gens

        self.assertEqual(group.order, 9)
        self.assertEqual((group.elt(first) * group.elt(second)).as_tuple, (1, 1))
        self.assertEqual((-group.elt(first)).as_tuple, (2, 0))
        self.assertEqual((group.elt(first) ** 3).as_tuple, group.one)

        heisenberg_group = heisenberg(3)
        for value in heisenberg_group.elts:
            self.assertEqual((heisenberg_group.elt(value) * (-heisenberg_group.elt(value))).as_tuple, heisenberg_group.one)

    def test_elementary_template(self):
        '''Check the elementary abelian cover template.'''
        template = elementary_template(5, 2)
        template_field, z, f, x, y = template.fct_field

        self.assertIsInstance(template, CoverTemplate)
        self.assertEqual(template.height, 2)
        self.assertEqual(template.group.order, 25)
        self.assertEqual(template.fcts, [template_field(f[0]), template_field(f[1])])
        self.assertEqual(template.gp_action[0], [template_field(z[0] + 1), template_field(z[1]), template_field(x), template_field(y)])
        self.assertEqual(template.gp_action[1], [template_field(z[0]), template_field(z[1] + 1), template_field(x), template_field(y)])

    def test_non_elementary_templates_smoke(self):
        '''Check representative non-elementary AS templates.'''
        templates = [
            witt_template(3, 2),
            heisenberg_template(3),
            quaternion_template(),
        ]

        self.assertEqual([(template.height, len(template.fcts), len(template.gp_action)) for template in templates], [
            (2, 2, 1),
            (3, 3, 3),
            (3, 3, 2),
        ])

    def test_elementary_cover_construction_and_alias(self):
        '''Check construction of an elementary AS cover.'''
        curve = self.quotient_curve()
        cover = elementary_cover([curve.x], prec=60)
        legacy_cover = as_cover(curve, elementary_template(5, 1), [curve.x], prec=60)

        self.assertIsInstance(cover, ArtinSchreierCover)
        self.assertIsInstance(legacy_cover, ArtinSchreierCover)
        self.assertEqual(cover.height, 1)
        self.assertEqual(cover.group.order, 5)
        self.assertEqual(cover.rhs[0], cover.x.function)
        self.assertEqual(len(cover.z), 1)

    def test_local_expansion_satisfies_as_equation(self):
        '''Check that z^p - z agrees with the substituted RHS locally.'''
        cover = self.elementary_cover()
        z_series = cover.z_series[0][0]
        x_series = cover.x_series[0]

        self.assertEqual(z_series ** cover.characteristic - z_series, x_series)
        self.assertEqual(cover.jumps[0][0], 2)

    def test_different_stabilizer_fiber_and_genus(self):
        '''Check basic ramification invariants of a one-generator cover.'''
        cover = self.elementary_cover()

        self.assertEqual(cover.exponent_of_different(0), 12)
        self.assertEqual(cover.stabilizer(0), cover.group.elts)
        self.assertEqual(cover.fiber(0), [cover.group.one])
        self.assertEqual(cover.genus(), 7)

    def test_function_arithmetic_reduction_and_group_action(self):
        '''Check function arithmetic, reduction, and generator action.'''
        cover = self.elementary_cover()
        generator = cover.group.gens[0]

        self.assertIsInstance(cover.z[0], ArtinSchreierFunction)
        self.assertEqual(cover.z[0].group_action(generator), cover.z[0] + 1)
        self.assertEqual(cover.z[0] ** cover.characteristic - cover.z[0], cover.x)
        self.assertEqual((cover.z[0] ** cover.characteristic).reduce(), cover.z[0] + cover.x)

    def test_forms_store_function_coefficients(self):
        '''Check that AS forms reuse AS functions as coefficients.'''
        cover = self.elementary_cover()
        omega = cover.z[0] * cover.dx

        self.assertIsInstance(omega, ArtinSchreierForm)
        self.assertIsInstance(omega.form, ArtinSchreierFunction)
        self.assertEqual(
            omega.expansion_at_infty(),
            cover.z[0].expansion_at_infty() * cover.dx_series[0],
        )

    def test_differential_and_cech_construction(self):
        '''Check differentiation and Cech-de Rham cocycle construction.'''
        cover = self.elementary_cover()
        differential = cover.z[0].diffn()
        cocycle = cover.cech(differential, cover.z[0])

        self.assertIsInstance(differential, ArtinSchreierForm)
        self.assertIsInstance(cocycle, ArtinSchreierDeRhamCocycle)
        self.assertEqual(cocycle.omega8, 0 * cover.dx)

    def test_group_action_matrices_and_polydifferentials(self):
        '''Check migrated group-action matrix and polydifferential helpers.'''
        cover = self.rational_elementary_cover()
        holomorphic_basis = cover.holomorphic_differentials_basis(threshold=5)
        holomorphic_matrices = cover.group_action_matrices_holo(basis=holomorphic_basis)
        poly_basis = cover.holo_polydifferentials_basis(2, threshold=5)
        poly_matrices = cover.group_action_matrices_poly(2, basis=poly_basis)

        self.assertEqual(len(holomorphic_matrices), 1)
        self.assertEqual(holomorphic_matrices[0].dimensions(), (1, 1))
        self.assertEqual(holomorphic_matrices[0][0, 0], 1)
        self.assertEqual(len(poly_basis), 1)
        self.assertIsInstance(poly_basis[0], ArtinSchreierPolyform)
        self.assertEqual(poly_matrices[0][0, 0], 1)
        self.assertIn("IndecomposableSummands", magma_module_decomposition(
            holomorphic_matrices[0],
            holomorphic_matrices[0],
            text=True,
        ))

    def test_riemann_roch_symmetric_and_dual_helpers(self):
        '''Check migrated Riemann-Roch, symmetric-product, and trace-dual helpers.'''
        cover = self.rational_elementary_cover()
        pole_orders = {(0, cover.group.one): 0}

        self.assertEqual(cover.riemann_roch_space(pole_orders, threshold=5), [cover.one])
        self.assertEqual(cover.riemann_roch_space_forms(pole_orders, threshold=5), [cover.dx])
        self.assertEqual(len(cover.symmetric_power_basis(2, threshold=5)), 1)
        self.assertEqual(cover.canonical_ideal_polynomials(2, threshold=5), [])

        dual = dual_elt(cover, cover.z[0])
        self.assertEqual((cover.z[0] * dual).trace(super=False), cover.one)
        self.assertEqual((cover.z[0] * dual.group_action(cover.group.gens[0])).trace(super=False), 0 * cover.x)
        omega = cover.holomorphic_differentials_basis(threshold=5)[0]
        self.assertEqual(combination_components(omega, cover.z[0], cover.z[0]), omega)


if __name__ == "__main__":
    unittest.main()
