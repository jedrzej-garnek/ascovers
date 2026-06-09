'''Cech-de Rham cocycles on Artin-Schreier type covers.'''

from itertools import product

from sage.all import prod, vector

from ascovers.as_covers.holomorphic import holomorphic_combinations_functions


class ArtinSchreierDeRhamCocycle:
    '''A Cech representative ``(omega0, f, omega8)`` of de Rham cohomology.'''

    def __init__(self, cover, omega, function):
        '''Create the cocycle with ``omega8 = omega - df``.'''
        self.curve = cover
        self.omega0 = omega if getattr(omega, "curve", None) is cover else cover.form(omega)
        self.f = function if getattr(function, "curve", None) is cover else cover.function(function)
        self.omega8 = self.omega0 - self.f.diffn()

    def __repr__(self):
        '''Return a compact textual representation.'''
        return f"( {self.omega0}, {self.f} )"

    def __eq__(self, other):
        '''Return whether two cocycles have the same components.'''
        return (
            isinstance(other, ArtinSchreierDeRhamCocycle)
            and self.curve is other.curve
            and self.omega0 == other.omega0
            and self.f == other.f
        )

    def __add__(self, other):
        '''Add two cocycles.'''
        if other.curve is not self.curve:
            raise TypeError("cannot add cocycles on different covers")
        return ArtinSchreierDeRhamCocycle(self.curve, self.omega0 + other.omega0, self.f + other.f)

    def __sub__(self, other):
        '''Subtract another cocycle.'''
        if other.curve is not self.curve:
            raise TypeError("cannot subtract cocycles on different covers")
        return ArtinSchreierDeRhamCocycle(self.curve, self.omega0 - other.omega0, self.f - other.f)

    def __neg__(self):
        '''Return the additive inverse.'''
        return ArtinSchreierDeRhamCocycle(self.curve, -self.omega0, -self.f)

    def __rmul__(self, scalar):
        '''Multiply the cocycle by a scalar.'''
        return ArtinSchreierDeRhamCocycle(self.curve, scalar * self.omega0, scalar * self.f)

    def reduce(self):
        '''Return a cocycle whose components have reduced representatives.'''
        return ArtinSchreierDeRhamCocycle(self.curve, self.omega0.reduce(), self.f.reduce())

    def group_action(self, element):
        '''Return the image of the cocycle under a group element.'''
        return ArtinSchreierDeRhamCocycle(
            self.curve,
            self.omega0.group_action(element),
            self.f.group_action(element),
        )

    def trace(self):
        '''Return the trace cocycle.'''
        return ArtinSchreierDeRhamCocycle(self.curve, self.omega0.trace(), self.f.trace())

    def coordinates(self, threshold=10, basis=0):
        '''Return coordinates in a de Rham basis when the required bases are known.'''
        cover = self.curve
        quotient_exponent = int(cover.quotient.exponent)
        quotient_degree = int(cover.quotient.polynomial.degree())
        function_field, polynomial_ring, _x_coordinate, _y_coordinate, _z_coordinates = cover.fct_field
        if basis == 0:
            basis = [
                cover.holomorphic_differentials_basis(),
                cover.cohomology_of_structure_sheaf_basis(),
                cover.de_rham_basis(threshold=threshold),
            ]
        holomorphic_forms, cohomology_basis, de_rham_basis = basis

        cohomology_coordinates = self.f.coordinates(basis=[holomorphic_forms, cohomology_basis])
        reduced = self
        for index in range(int(cover.genus())):
            reduced -= cohomology_coordinates[index] * de_rham_basis[index + int(cover.genus())]
        coordinate_prefix = [cover.base_ring.zero()] * int(cover.genus()) + list(cohomology_coordinates)

        if reduced.f.function in polynomial_ring:
            reduced.omega0 -= reduced.f.diffn()
            return vector(cover.base_ring, coordinate_prefix) + vector(
                cover.base_ring,
                list(reduced.omega0.coordinates(basis=holomorphic_forms))
                + [cover.base_ring.zero()] * int(cover.genus()),
            )

        candidate_items = []
        z_exponent_ranges = [range(int(cover.characteristic)) for _index in range(cover.height)]
        for x_power in range(0, int(threshold) * quotient_degree):
            for y_power in range(0, quotient_exponent):
                for z_exponents in product(*z_exponent_ranges):
                    expression = cover.x ** x_power
                    expression *= prod(
                        cover.z[index] ** int(z_exponents[index])
                        for index in range(cover.height)
                    )
                    expression *= cover.y ** y_power
                    candidate_items.append((expression, expression.expansion_at_infty()))
        candidate_items.append((reduced.f, reduced.f.expansion_at_infty()))
        candidates = holomorphic_combinations_functions(candidate_items, 0)

        for place in range(int(cover.nb_of_pts_at_infty)):
            for element in cover.fiber(place=place):
                if place != 0 or element != cover.group.one:
                    candidates = [
                        (function, function.group_action(element).expansion_at_infty(place=place))
                        for function in candidates
                    ]
                    candidates = holomorphic_combinations_functions(candidates, 0)

        for function in candidates:
            if function.function not in polynomial_ring:
                for scalar in cover.base_ring:
                    if reduced.f.function - scalar * function.function in polynomial_ring:
                        reduced.f = cover.function(reduced.f.function - scalar * function.function)
                        return vector(cover.base_ring, coordinate_prefix) + vector(
                            cover.base_ring,
                            reduced.coordinates(threshold=threshold, basis=basis),
                        )
        raise ValueError("increase threshold")

    def pairing(self, other):
        '''Return the de Rham cup-product pairing with another cocycle.'''
        cover = self.curve
        omega1 = self.omega0
        function1 = self.f
        omega2 = other.omega0
        function2 = other.f
        return (
            omega1.serre_duality_pairing(function2)
            - omega2.serre_duality_pairing(function1)
            - sum(
                (function2 * function1.diffn()).residue(place=place)
                for place in range(int(cover.nb_of_pts_at_infty))
            )
        )


ArtinSchreierCechCocycle = ArtinSchreierDeRhamCocycle
ArtinSchreierCech = ArtinSchreierDeRhamCocycle
as_cech = ArtinSchreierDeRhamCocycle
