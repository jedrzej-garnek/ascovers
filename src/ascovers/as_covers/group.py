'''Finite groups used as Galois groups of Artin-Schreier type covers.'''

from itertools import product

from sage.all import ZZ


class FiniteGroup:
    '''A small finite group described by explicit elements and operations.'''

    def __init__(self, name, short_name, elements, identity, multiply, inverse, generators):
        '''Create a finite group from explicit multiplication and inversion maps.

        INPUT:

        - ``name`` -- human-readable group name.
        - ``short_name`` -- compact name used in cover representations.
        - ``elements`` -- iterable of raw element labels.
        - ``identity`` -- raw label of the identity element.
        - ``multiply`` -- function multiplying two raw labels.
        - ``inverse`` -- function inverting a raw label.
        - ``generators`` -- raw labels of generators, in template-action order.
        '''
        self.name = name
        self.short_name = short_name
        self.elts = list(elements)
        self.one = identity
        self.mult = multiply
        self.inv = inverse
        self.gens = list(generators)
        self.order = len(self.elts)

    def __repr__(self):
        '''Return the group name.'''
        return self.name

    def element(self, value):
        '''Return the group element represented by ``value``.'''
        return FiniteGroupElement(value, self)

    def elt(self, value):
        '''Compatibility alias for ``element``.'''
        return self.element(value)

    def one_element(self):
        '''Return the identity as a group-element object.'''
        return self.element(self.one)

    def ONE(self):
        '''Compatibility alias returning the identity element.'''
        return self.one_element()

    def generators(self):
        '''Return the distinguished generators as group-element objects.'''
        return [self.element(value) for value in self.gens]

    def GENS(self):
        '''Compatibility alias returning the distinguished generators.'''
        return self.generators()


class FiniteGroupElement:
    '''An element of an explicitly described finite group.'''

    def __init__(self, value, group):
        '''Create the element represented by ``value`` in ``group``.'''
        self.as_tuple = value
        self.group = group

    def __repr__(self):
        '''Return the raw label representation.'''
        return str(self.as_tuple)

    def __eq__(self, other):
        '''Return whether two group elements have the same parent and label.'''
        return (
            isinstance(other, FiniteGroupElement)
            and self.group is other.group
            and self.as_tuple == other.as_tuple
        )

    def __hash__(self):
        '''Return a hash compatible with equality.'''
        return hash((id(self.group), self.as_tuple))

    def __mul__(self, other):
        '''Multiply two group elements.'''
        if isinstance(other, FiniteGroupElement):
            if other.group is not self.group:
                raise TypeError("cannot multiply elements from different groups")
            other = other.as_tuple
        return FiniteGroupElement(self.group.mult(self.as_tuple, other), self.group)

    def __rmul__(self, other):
        '''Multiply a raw label by this element.'''
        return FiniteGroupElement(self.group.mult(other, self.as_tuple), self.group)

    def __neg__(self):
        '''Return the inverse element.

        The legacy Sage code used unary ``-`` for inversion, and the alias is kept
        for compatibility with old notebooks.
        '''
        return FiniteGroupElement(self.group.inv(self.as_tuple), self.group)

    def inverse(self):
        '''Return the inverse element.'''
        return -self

    def __pow__(self, exponent):
        '''Return this element raised to an integer power.'''
        exponent = ZZ(exponent)
        if exponent == 0:
            return self.group.ONE()
        if exponent < 0:
            return (-self) ** (-exponent)

        result = self.group.ONE()
        base = self
        while exponent:
            if exponent % 2:
                result = result * base
            base = base * base
            exponent //= 2
        return result


def cyclic_group(characteristic, exponent):
    '''Return the cyclic group of order ``characteristic^exponent``.'''
    order = int(characteristic) ** int(exponent)
    return FiniteGroup(
        name=f"cyclic group of order {characteristic}^{exponent}",
        short_name="Z/p^n",
        elements=range(order),
        identity=0,
        multiply=lambda first, second: (first + second) % order,
        inverse=lambda value: (-value) % order,
        generators=[1],
    )


def elementary_abelian_group(characteristic, rank):
    '''Return the elementary abelian group ``(Z/p)^rank``.'''
    characteristic = int(characteristic)
    rank = int(rank)
    elements = [tuple(value) for value in product(range(characteristic), repeat=rank)]
    generators = []
    for index in range(rank):
        generator = [0] * rank
        generator[index] = 1
        generators.append(tuple(generator))

    return FiniteGroup(
        name=f"(Z/{characteristic})^{rank}",
        short_name=f"(Z/{characteristic})^{rank}",
        elements=elements,
        identity=tuple([0] * rank),
        multiply=lambda first, second: tuple(
            (first[index] + second[index]) % characteristic for index in range(rank)
        ),
        inverse=lambda value: tuple((-value[index]) % characteristic for index in range(rank)),
        generators=generators,
    )


def heisenberg_group(characteristic):
    '''Return the exponent-p Heisenberg group of order ``p^3``.'''
    characteristic = int(characteristic)
    elements = [
        (first, second, third)
        for first in range(characteristic)
        for second in range(characteristic)
        for third in range(characteristic)
    ]

    def multiply(first, second):
        return (
            (first[0] + second[0]) % characteristic,
            (first[1] + second[1]) % characteristic,
            (-first[1] * second[0] + first[2] + second[2]) % characteristic,
        )

    def inverse(value):
        return (
            (-value[0]) % characteristic,
            (-value[1]) % characteristic,
            (-value[2] - value[0] * value[1]) % characteristic,
        )

    return FiniteGroup(
        name=f"Heisenberg group E({characteristic}^3)",
        short_name="E(p^3)",
        elements=elements,
        identity=(0, 0, 0),
        multiply=multiply,
        inverse=inverse,
        generators=[(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    )


def quaternion_multiply(first, second):
    '''Multiply raw labels in the legacy model of Q8.'''
    result = [
        (first[0] + second[0] + 2 * first[1] * second[0]) % 4,
        (first[1] + second[1]) % 4,
    ]
    if result[1] in (2, 3):
        result[0] = (result[0] + 2) % 4
        result[1] = (result[1] - 2) % 4
    return tuple(result)


def quaternion_inverse(value):
    '''Invert a raw label in the legacy model of Q8.'''
    result = [
        ((-1) ** (value[0] * value[1]) * (-value[0])) % 4,
        (-value[1]) % 4,
    ]
    if result[1] in (2, 3):
        result[0] = (result[0] + 2) % 4
        result[1] = (result[1] - 2) % 4
    return tuple(result)


def quaternion_group():
    '''Return the quaternion group Q8.'''
    return FiniteGroup(
        name="Q8",
        short_name="Q8",
        elements=[(first, second) for first in range(4) for second in range(2)],
        identity=(0, 0),
        multiply=quaternion_multiply,
        inverse=quaternion_inverse,
        generators=[(1, 0), (0, 1)],
    )


def hypoelementary_multiply(characteristic, order, character_value, first_a, first_b, second_a, second_b):
    '''Multiply raw labels in a hypoelementary group.'''
    return ((first_a + second_a) % order, (character_value ** second_a * first_b + second_b) % characteristic)


def hypoelementary_inverse(characteristic, order, character_value, first_a, first_b):
    '''Invert a raw label in a hypoelementary group.'''
    return hypoelementary_multiply(characteristic, order, character_value, 0, characteristic - first_b, order - first_a, 0)


def hypoelementary_group(characteristic, order, character_value):
    '''Return ``Z/p semidirect Z/m`` for a specified character value.'''
    return FiniteGroup(
        name=(
            f"Hypoelementary group Z/{characteristic} semidirect Z/{order}, "
            f"with character value {character_value}"
        ),
        short_name=f"Z/{characteristic} semidirect Z/{order}",
        elements=[(first, second) for first in range(int(order)) for second in range(int(characteristic))],
        identity=(0, 0),
        multiply=lambda first, second: hypoelementary_multiply(
            characteristic, order, character_value, first[0], first[1], second[0], second[1]
        ),
        inverse=lambda value: hypoelementary_inverse(characteristic, order, character_value, value[0], value[1]),
        generators=[(1, 0), (0, 1)],
    )


def dihedral_group(characteristic, exponent):
    '''Return the dihedral group of order ``2 * characteristic^exponent``.'''
    modulus = int(characteristic) ** int(exponent)
    return FiniteGroup(
        name=f"Dihedral group D({modulus})",
        short_name="D(p^n)",
        elements=[(first, second) for first in range(modulus) for second in range(2)],
        identity=(0, 0),
        multiply=lambda first, second: (
            (first[0] + (-1) ** first[1] * second[0]) % modulus,
            (first[1] + second[1]) % 2,
        ),
        inverse=lambda value: (((-1) ** value[1] * (-value[0])) % modulus, value[1]),
        generators=[(1, 0), (0, 1)],
    )

def A4_group():
    '''Return the alternating group A_4.

    Presentation:
        <sigma, rho | sigma^2 = rho^3 = (sigma rho)^3 = 1>

    Elements are represented as triples (i, j, k), corresponding to
        sigma^i rho^j (sigma rho)^k.
    '''
    def compose(first, second):
        # Permutations are tuples, with value[i] = image of i.
        # This returns first ∘ second.
        return tuple(first[i] for i in second)

    def invert_perm(value):
        inverse_value = [0] * len(value)
        for index, image in enumerate(value):
            inverse_value[image] = index
        return tuple(inverse_value)

    identity_perm = (0, 1, 2, 3)

    sigma_perm = (1, 0, 3, 2)  # (1 2)(3 4)
    rho_perm = (1, 2, 0, 3)    # (1 2 3)

    sigma_rho_perm = compose(sigma_perm, rho_perm)

    def power(value, exponent):
        result = identity_perm
        for _ in range(exponent):
            result = compose(result, value)
        return result

    def to_perm(value):
        i, j, k = value
        return compose(
            compose(power(sigma_perm, i), power(rho_perm, j)),
            power(sigma_rho_perm, k),
        )

    elements = [
        (i, j, k)
        for i in range(2)
        for j in range(3)
        for k in range(2)
    ]

    perm_to_element = {
        to_perm(element): element
        for element in elements
    }

    return FiniteGroup(
        name="Alternating group A_4",
        short_name="A_4",
        elements=elements,
        identity=(0, 0, 0),
        multiply=lambda first, second: perm_to_element[
            compose(to_perm(first), to_perm(second))
        ],
        inverse=lambda value: perm_to_element[
            invert_perm(to_perm(value))
        ],
        generators=[
            (1, 0, 0),  # sigma
            (0, 1, 0),  # rho
        ],
    )

group = FiniteGroup
group_elt = FiniteGroupElement
cyclic_gp = cyclic_group
elementary_gp = elementary_abelian_group
heisenberg = heisenberg_group
quaternion_gp = quaternion_group
hypoelementary = hypoelementary_group
dihedreal_gp = dihedral_group