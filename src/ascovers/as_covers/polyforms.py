'''Polydifferentials and canonical-ideal helpers for AS covers.'''

from itertools import product

from sage.all import LCM, PolynomialRing, ZZ, matrix, prod, vector

from ascovers.as_covers.reduction import artin_schreier_reduction
from ascovers.linear_algebra import (
    linear_combination_coefficients,
    linear_representation_polynomials,
)


class ArtinSchreierPolyform:
    '''An element of ``H^0(Omega^tensor mult)`` represented as ``h dx^mult``.'''

    def __init__(self, coefficient, mult):
        '''Create a polydifferential from an AS function coefficient.'''
        self.form = coefficient
        self.curve = coefficient.curve
        self.mult = int(mult)

    def __repr__(self):
        '''Return a textual representation of the polydifferential.'''
        return f"({self.form}) dx^tensor {self.mult}"

    def __add__(self, other):
        '''Add two polydifferentials of the same tensor degree.'''
        if other.curve is not self.curve or other.mult != self.mult:
            raise TypeError("cannot add polydifferentials on different spaces")
        return ArtinSchreierPolyform(self.form + other.form, self.mult)

    def __sub__(self, other):
        '''Subtract another polydifferential.'''
        if other.curve is not self.curve or other.mult != self.mult:
            raise TypeError("cannot subtract polydifferentials on different spaces")
        return ArtinSchreierPolyform(self.form - other.form, self.mult)

    def __neg__(self):
        '''Return the additive inverse.'''
        return ArtinSchreierPolyform(-self.form, self.mult)

    def __rmul__(self, scalar):
        '''Multiply by a scalar.'''
        return ArtinSchreierPolyform(scalar * self.form, self.mult)

    def expansion_at_infty(self, place=0, prec=None):
        '''Return the Laurent expansion at a place above infinity.'''
        return self.form.expansion_at_infty(place=place, prec=prec) * (
            self.curve.dx.expansion_at_infty(place=place, prec=prec) ** self.mult
        )

    def reduce(self):
        '''Return this polydifferential with reduced coefficient.'''
        return ArtinSchreierPolyform(self.form.reduce(), self.mult)

    def coordinates(self, basis=0):
        '''Return coordinates in a basis of holomorphic polydifferentials.'''
        if basis == 0:
            basis = self.curve.holo_polydifferentials_basis(mult=self.mult)
        _function_field, polynomial_ring, _x_coordinate, _y_coordinate, _z_coordinates = self.curve.fct_field
        common_denominator = LCM([element.form.function.denominator() for element in basis])
        basis_polynomials = [common_denominator * element.form.function for element in basis]
        reduced = self.reduce()
        target_polynomial = common_denominator * reduced.form.function
        return linear_representation_polynomials(
            polynomial_ring(target_polynomial),
            [polynomial_ring(polynomial) for polynomial in basis_polynomials],
        )

    def group_action(self, element):
        '''Return the image under a group element.'''
        return ArtinSchreierPolyform(self.form.group_action(element), self.mult)


def holomorphic_polydifferentials_basis(cover, mult, threshold=8):
    '''Return a searched basis of ``H^0(Omega^tensor mult)``.'''
    tensor_power = int(mult)
    dx_valuation = cover.dx.valuation()
    functions = cover.at_most_poles(tensor_power * dx_valuation, threshold=threshold)
    result = [ArtinSchreierPolyform(function, tensor_power) for function in functions]
    if tensor_power == 1 and len(result) < cover.genus():
        raise ValueError("increase threshold; not all polydifferentials were found")
    if cover.genus() > 1 and tensor_power > 1:
        expected_dimension = (2 * tensor_power - 1) * (cover.genus() - 1)
        if len(result) < expected_dimension:
            raise ValueError("increase threshold; not all polydifferentials were found")
    return result


def _is_non_decreasing(values):
    '''Return whether a tuple is nondecreasing.'''
    return all(first <= second for first, second in zip(values, values[1:]))


def symmetric_power_basis(cover, tensor_power, threshold=8):
    '''Return the natural monomial basis of ``Sym^tensor_power H^0(Omega)``.'''
    genus = int(cover.genus())
    holomorphic_basis = cover.holomorphic_differentials_basis(threshold=threshold)
    index_ranges = [range(genus) for _index in range(int(tensor_power))]
    result = []
    for indices in product(*index_ranges):
        if _is_non_decreasing(indices):
            result.append(
                SymmetricProductOfForms([[1] + [holomorphic_basis[index] for index in indices]])
            )
    return result


def canonical_ideal(cover, tensor_power, threshold=8):
    '''Return a basis for the degree-``tensor_power`` canonical ideal component.'''
    base_field = cover.base_ring
    symmetric_basis = symmetric_power_basis(cover, tensor_power, threshold=threshold)
    poly_basis = cover.holo_polydifferentials_basis(tensor_power, threshold=threshold)
    if not symmetric_basis:
        return []
    if not poly_basis:
        return symmetric_basis

    multiplication_matrix = matrix(base_field, len(poly_basis), len(symmetric_basis))
    for column, symmetric_product in enumerate(symmetric_basis):
        multiplication_matrix[:, column] = vector(
            base_field,
            symmetric_product.multiply().coordinates(basis=poly_basis),
        )

    result = []
    for relation in multiplication_matrix.right_kernel().basis():
        relation_product = 0 * symmetric_basis[0]
        for index, coefficient in enumerate(relation):
            relation_product += coefficient * symmetric_basis[index]
        result.append(relation_product)
    return result


def canonical_ideal_polynomials(cover, tensor_power, threshold=8):
    '''Return polynomial equations for a canonical ideal component.'''
    return [relation.polynomial() for relation in canonical_ideal(cover, tensor_power, threshold=threshold)]


class SymmetricProductOfForms:
    '''A linear combination of symmetric products of holomorphic forms.'''

    def __init__(self, forms_and_coefficients):
        '''Create from rows ``[coefficient, form_1, ..., form_n]``.'''
        simplified = []
        for entry in forms_and_coefficients:
            coefficient = entry[0]
            factors = entry[1:]
            existing_factors = [item[1:] for item in simplified]
            if coefficient != 0 and factors not in existing_factors:
                simplified.append(list(entry))
            elif factors in existing_factors:
                index = existing_factors.index(factors)
                simplified[index][0] += coefficient

        if not simplified:
            simplified = [[0] + forms_and_coefficients[0][1:]]

        self.tuples = simplified
        self.n = len(simplified[0]) - 1
        self.curve = simplified[0][1].curve

    def __repr__(self):
        '''Return a textual representation of this symmetric product.'''
        terms = []
        for entry in self.tuples:
            if entry[0] == 0:
                continue
            factors = " * ".join(f"({factor})" for factor in entry[1:])
            terms.append(f"{entry[0]} * {factors}")
        return " + ".join(terms) if terms else "0"

    def __add__(self, other):
        '''Add two symmetric-product combinations.'''
        return SymmetricProductOfForms(self.tuples + other.tuples)

    def __sub__(self, other):
        '''Subtract another symmetric-product combination.'''
        return self + (-1) * other

    def __neg__(self):
        '''Return the additive inverse.'''
        return (-1) * self

    def __rmul__(self, scalar):
        '''Multiply by a scalar.'''
        if isinstance(scalar, int) or scalar in ZZ or scalar in self.curve.base_ring:
            return SymmetricProductOfForms(
                [[scalar * entry[0]] + entry[1:] for entry in self.tuples]
            )
        return NotImplemented

    def coordinates(self, basis=0, dictionary=False):
        '''Return tensor coordinates in the holomorphic-form basis.'''
        cover = self.curve
        genus = int(cover.genus())
        if basis == 0:
            basis = cover.holomorphic_differentials_basis()
        index_ranges = [range(genus) for _index in range(self.n)]
        result = {indices: cover.base_ring.zero() for indices in product(*index_ranges)}

        for entry in self.tuples:
            coordinate_lists = [form.coordinates(basis=basis) for form in entry[1:]]
            for indices in product(*index_ranges):
                coefficient = cover.base_ring.one()
                for factor_index in range(self.n):
                    coefficient *= coordinate_lists[factor_index][indices[factor_index]]
                result[indices] += entry[0] * coefficient

        if dictionary:
            return result
        return [result[indices] for indices in product(*index_ranges)]

    def canonical_ideal_coordinates(self, basis=0):
        '''Return coordinates in a canonical-ideal basis.'''
        cover = self.curve
        if basis == 0:
            holomorphic_basis = cover.holomorphic_differentials_basis()
            ideal_basis = [
                relation.coordinates(basis=holomorphic_basis)
                for relation in cover.canonical_ideal(2)
            ]
        else:
            holomorphic_basis, ideal_basis = basis
        return linear_combination_coefficients(
            self.coordinates(basis=holomorphic_basis),
            ideal_basis,
        )

    def multiply(self):
        '''Apply ``Sym^n H^0(Omega) -> H^0(Omega^tensor n)``.'''
        cover = self.curve
        result = ArtinSchreierPolyform(0 * cover.x, self.n)
        for entry in self.tuples:
            coefficient = cover.one
            for differential in entry[1:]:
                coefficient *= differential.form
            reduced = cover.function(
                artin_schreier_reduction(cover, coefficient.function)
            )
            result += ArtinSchreierPolyform(entry[0] * reduced, self.n)
        return result

    def polynomial(self):
        '''Return the corresponding polynomial in canonical coordinates.'''
        cover = self.curve
        genus = int(cover.genus())
        coordinate_dictionary = self.coordinates(dictionary=True)
        polynomial_ring = PolynomialRing(cover.base_ring, "X", genus)
        variables = polynomial_ring.gens()
        result = polynomial_ring.zero()
        index_ranges = [range(genus) for _index in range(self.n)]
        for indices in product(*index_ranges):
            monomial = polynomial_ring.one()
            for factor_index in range(self.n):
                monomial *= variables[indices[factor_index]]
            result += coordinate_dictionary[indices] * monomial
        return result

    def group_action(self, element):
        '''Return the image under a group element.'''
        return SymmetricProductOfForms(
            [
                [entry[0]] + [factor.group_action(element) for factor in entry[1:]]
                for entry in self.tuples
            ]
        )


def group_action_canonical_ideal(cover, tensor_power, threshold=8):
    '''Return group-action matrices on a canonical-ideal component.'''
    ideal_basis = canonical_ideal(cover, tensor_power, threshold=threshold)
    if not ideal_basis:
        return []
    base_field = cover.base_ring
    basis_polynomials = [relation.polynomial() for relation in ideal_basis]
    result = []
    for group_element in cover.group.gens:
        action_matrix = matrix(base_field, len(ideal_basis), len(ideal_basis))
        for column, relation in enumerate(ideal_basis):
            transformed = relation.group_action(group_element).polynomial()
            action_matrix[:, column] = vector(
                base_field,
                linear_representation_polynomials(transformed, basis_polynomials),
            )
        result.append(action_matrix)
    return result


as_polyform = ArtinSchreierPolyform
as_holo_polydifferentials_basis = holomorphic_polydifferentials_basis
as_symmetric_product_forms = SymmetricProductOfForms
as_symmetric_power_basis = symmetric_power_basis
as_canonical_ideal = canonical_ideal
as_canonical_ideal_polynomials = canonical_ideal_polynomials
as_matrices_group_action_canonical_ideal = group_action_canonical_ideal
