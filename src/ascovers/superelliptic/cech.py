'''Cech-de Rham cocycles on superelliptic curves.'''

from sage.all import vector

from ascovers.superelliptic.decomposition import decomposition_g0_g8
from ascovers.superelliptic.form import SuperellipticForm
from ascovers.superelliptic.function import SuperellipticFunction


class SuperellipticDeRhamCocycle:
    '''A Cech-de Rham cocycle representing a class in algebraic de Rham cohomology.'''

    def __init__(self, curve, omega0, transition_function):
        '''Create the cocycle represented by (omega0, f, omega_infinity).

        The cocycle is stored using the convention

        ``omega_infinity = omega0 - d(transition_function)``.
        '''
        self.curve = curve
        self.omega0 = self._coerce_form(omega0)
        self.f = self._coerce_function(transition_function)
        self.omega8 = self.omega0 - self.f.diffn()
        self.transition_function = self.f
        self.omega_infinity = self.omega8

    def _coerce_form(self, value):
        '''Coerce a value to a differential form on this cocycle's curve.'''
        if isinstance(value, SuperellipticForm):
            if value.curve is not self.curve:
                raise TypeError("cannot use a form from a different curve")
            return value
        return SuperellipticForm(self.curve, value)

    def _coerce_function(self, value):
        '''Coerce a value to a function on this cocycle's curve.'''
        if isinstance(value, SuperellipticFunction):
            if value.curve is not self.curve:
                raise TypeError("cannot use a function from a different curve")
            return value
        return SuperellipticFunction(self.curve, value)

    def __eq__(self, other):
        '''Return whether two cocycles have the same stored representatives.'''
        if not isinstance(other, SuperellipticDeRhamCocycle):
            return False
        return self.curve is other.curve and self.omega0 == other.omega0 and self.f == other.f

    def __add__(self, other):
        '''Add two Cech-de Rham cocycles on the same curve.'''
        self._check_same_curve(other)
        return SuperellipticDeRhamCocycle(self.curve, self.omega0 + other.omega0, self.f + other.f)

    def __sub__(self, other):
        '''Subtract another Cech-de Rham cocycle on the same curve.'''
        self._check_same_curve(other)
        return SuperellipticDeRhamCocycle(self.curve, self.omega0 - other.omega0, self.f - other.f)

    def __neg__(self):
        '''Return the additive inverse of this Cech-de Rham cocycle.'''
        return SuperellipticDeRhamCocycle(self.curve, -self.omega0, -self.f)

    def __rmul__(self, scalar):
        '''Multiply this cocycle by a scalar.'''
        return SuperellipticDeRhamCocycle(self.curve, scalar * self.omega0, scalar * self.f)

    def __repr__(self):
        '''Return a readable representation of the Cech-de Rham triple.'''
        return "(" + str(self.omega0) + ", " + str(self.f) + ", " + str(self.omega8) + ")"

    def _check_same_curve(self, other):
        '''Raise an error if two cocycles live on different curves.'''
        if not isinstance(other, SuperellipticDeRhamCocycle):
            raise TypeError("expected a SuperellipticDeRhamCocycle")
        if other.curve is not self.curve:
            raise TypeError("cannot combine cocycles on different curves")

    def verschiebung(self):
        '''Return the Verschiebung image of this cocycle.'''
        return SuperellipticDeRhamCocycle(
            self.curve,
            self.omega0.cartier(),
            self.curve.function(self.curve.function_field.zero()),
        )

    def frobenius(self):
        '''Return the Frobenius image of this cocycle.'''
        characteristic = self.curve.characteristic
        if characteristic == 0:
            raise ValueError("Frobenius is only implemented in positive characteristic")
        return SuperellipticDeRhamCocycle(
            self.curve,
            self.curve.form(self.curve.function_field.zero()),
            self.f ** characteristic,
        )

    def coordinates(self, basis=0):
        '''Return coordinates in the standard de Rham basis.'''
        genus = int(self.curve.genus())
        base_field = self.curve.base_ring
        if basis == 0:
            basis_holomorphic = self.curve.holomorphic_differentials_basis()
            basis_structure_sheaf = self.curve.cohomology_of_structure_sheaf_basis()
            basis_de_rham = self.curve.de_rham_basis()
        else:
            basis_holomorphic, basis_structure_sheaf, basis_de_rham = basis

        for index, basis_element in enumerate(basis_de_rham):
            if self == basis_element:
                coordinates = [base_field.zero()] * (2 * genus)
                coordinates[index] = base_field.one()
                return vector(base_field, coordinates)

        zero_function = self.curve.function(self.curve.function_field.zero())
        zero_form = self.curve.form(self.curve.function_field.zero())

        if self.f == zero_function and self.omega0 == zero_form:
            return vector(base_field, [base_field.zero()] * (2 * genus))

        if self.f == zero_function:
            coordinates = list(self.omega0.coordinates(basis=basis_holomorphic))
            coordinates += [base_field.zero()] * genus
            return vector(base_field, coordinates)

        if self.omega0 == zero_form:
            coordinates = [base_field.zero()] * genus
            coordinates += list(self.f.coordinates(basis=basis_holomorphic))
            return vector(base_field, coordinates)

        structure_coordinates = self.f.coordinates(basis=basis_holomorphic)
        reduced = self
        for index in range(genus):
            reduced -= structure_coordinates[index] * basis_de_rham[index + genus]

        affine_part, _infinity_part, _cohomology_part = decomposition_g0_g8(reduced.f)
        holomorphic_cocycle = SuperellipticDeRhamCocycle(
            self.curve,
            reduced.omega0 - affine_part.diffn(),
            zero_function,
        )
        coordinates = [base_field.zero()] * genus + list(structure_coordinates)
        return vector(base_field, coordinates) + holomorphic_cocycle.coordinates(basis=basis)

    def is_cocycle(self):
        '''Return True if the stored triple satisfies the regularity checks.'''
        omega0_regular = self.omega0.is_regular_on_U0()
        omega_infinity_regular = self.omega8.is_regular_on_Uinfty()

        if omega0_regular and omega_infinity_regular:
            return True
        if not omega0_regular and not omega_infinity_regular:
            return "omega0 & omega_infinity"
        if not omega0_regular:
            return "omega0"
        return "omega_infinity"

    def is_valid(self):
        '''Return whether this object is a valid Cech-de Rham cocycle.'''
        return self.is_cocycle() is True


SuperellipticCechCocycle = SuperellipticDeRhamCocycle
SuperellipticCech = SuperellipticDeRhamCocycle
superelliptic_cech = SuperellipticDeRhamCocycle
