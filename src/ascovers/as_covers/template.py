'''Equation templates for Artin-Schreier type p-group covers.'''

from sage.all import GF, PolynomialRing, QQ

from ascovers.as_covers.group import (
    cyclic_group,
    dihedral_group,
    elementary_abelian_group,
    heisenberg_group,
    hypoelementary_group,
    quaternion_group,
    A4_group,
)


class CoverTemplate:
    '''A template for a p-group cover over a superelliptic quotient curve.

    A template stores symbolic right-hand sides for equations in generators
    ``z_i`` and symbolic group actions on the variables ``z_i, x, y``.  Concrete
    covers are obtained by substituting functions on the quotient curve for the
    formal parameters ``f_i``.
    '''

    def __init__(self, height, field, group, equations, generator_actions):
        '''Create a cover template.

        ``equations`` and ``generator_actions`` may either be lists of symbolic
        expressions or callables accepting ``(z, f, x, y)`` and returning those
        lists in the template polynomial ring.
        '''
        self.height = int(height)
        self.field = field
        self.group = group

        variable_names = [f"z{index}" for index in range(self.height)]
        variable_names += [f"f{index}" for index in range(self.height)]
        variable_names += ["x", "y"]
        polynomial_ring = PolynomialRing(field, tuple(variable_names))
        generators = polynomial_ring.gens()
        z_generators = generators[: self.height]
        f_generators = generators[self.height : 2 * self.height]
        x_generator = generators[-2]
        y_generator = generators[-1]
        function_field = polynomial_ring.fraction_field()

        if callable(equations):
            equations = equations(z_generators, f_generators, x_generator, y_generator)
        if callable(generator_actions):
            generator_actions = generator_actions(z_generators, f_generators, x_generator, y_generator)

        self.fct_field = (
            function_field,
            z_generators,
            f_generators,
            x_generator,
            y_generator,
        )
        self.fcts = [function_field(equation) for equation in equations]
        self.gp_action = [
            [function_field(entry) for entry in action]
            for action in generator_actions
        ]

        if len(self.fcts) != self.height:
            raise ValueError("a template must have one equation for each AS generator")
        for action in self.gp_action:
            if len(action) != self.height + 2:
                raise ValueError("each generator action must specify z_i, x, and y")


def elementary_template(characteristic, rank):
    '''Return the template for an elementary abelian ``(Z/p)^rank`` cover.'''
    cover_group = elementary_abelian_group(characteristic, rank)

    def equations(_z, f, _x, _y):
        return [f[index] for index in range(int(rank))]

    def actions(z, _f, x, y):
        return [
            [z[index] + (generator_index == index) for index in range(int(rank))] + [x, y]
            for generator_index in range(int(rank))
        ]

    return CoverTemplate(rank, GF(characteristic), cover_group, equations, actions)


def elementary_cover(functions, prec=10):
    '''Return the elementary abelian cover defined by the given quotient functions.'''
    from ascovers.as_covers.cover import ArtinSchreierCover

    quotient = functions[0].curve
    return ArtinSchreierCover(
        quotient,
        elementary_template(quotient.characteristic, len(functions)),
        functions,
        branch_points=[],
        prec=prec,
    )


def heisenberg_template(characteristic):
    '''Return the template for the exponent-p Heisenberg cover.'''
    cover_group = heisenberg_group(characteristic)

    def equations(z, f, _x, _y):
        return [f[0], f[1], f[2] + (z[0] - z[1]) * f[1]]

    def actions(z, _f, x, y):
        return [
            [z[0] + 1, z[1], z[2] + z[1], x, y],
            [z[0] + 1, z[1] + 1, z[2], x, y],
            [z[0], z[1], z[2] - 1, x, y],
        ]

    return CoverTemplate(3, GF(characteristic), cover_group, equations, actions)


def heisenberg_cover(functions, prec=10):
    '''Return the Heisenberg cover defined by three quotient functions.'''
    from ascovers.as_covers.cover import ArtinSchreierCover

    quotient = functions[0].curve
    return ArtinSchreierCover(
        quotient,
        heisenberg_template(quotient.characteristic),
        functions,
        branch_points=[],
        prec=prec,
    )


def witt_polynomial(variables, characteristic, index):
    '''Return the Witt ghost polynomial used by the legacy template formulas.'''
    return sum(
        characteristic ** degree * variables[degree] ** (characteristic ** (index - degree - 1))
        for degree in range(index)
    )


def witt_sum(characteristic, index):
    '''Return the universal Witt addition polynomial in coordinate ``index``.'''
    variable_names = [f"X{degree}" for degree in range(index + 1)]
    variable_names += [f"Y{degree}" for degree in range(index + 1)]
    polynomial_ring = PolynomialRing(QQ, tuple(variable_names))
    variables = polynomial_ring.gens()
    x_variables = variables[: index + 1]
    y_variables = variables[index + 1 :]

    if index == 0:
        return x_variables[0] + y_variables[0]

    previous_sums = []
    for previous_index in range(index):
        previous = witt_sum(characteristic, previous_index)
        previous_ring = previous.parent()
        previous_x = previous_ring.gens()[: previous_index + 1]
        previous_y = previous_ring.gens()[previous_index + 1 :]
        substitutions = {
            previous_x[degree]: x_variables[degree]
            for degree in range(previous_index + 1)
        }
        substitutions.update(
            {
                previous_y[degree]: y_variables[degree]
                for degree in range(previous_index + 1)
            }
        )
        previous_sums.append(previous.substitute(substitutions))

    numerator = witt_polynomial(x_variables, characteristic, index + 1)
    numerator += witt_polynomial(y_variables, characteristic, index + 1)
    numerator -= sum(
        characteristic ** degree
        * previous_sums[degree] ** (characteristic ** (index - degree))
        for degree in range(index)
    )
    return QQ(1) / (characteristic ** index) * numerator


def witt_sum_mod_p(characteristic, index):
    '''Return the Witt addition polynomial in characteristic ``p``.'''
    variable_names = [f"X{degree}" for degree in range(index + 1)]
    variable_names += [f"Y{degree}" for degree in range(index + 1)]
    rational_polynomial = witt_sum(characteristic, index)
    polynomial_ring = PolynomialRing(GF(characteristic), tuple(variable_names))
    return polynomial_ring(rational_polynomial)


def witt_template(characteristic, height):
    '''Return the Artin-Schreier-Witt template of length ``height``.'''
    cover_group = cyclic_group(characteristic, height)

    def equations(z, f, _x, _y):
        rhs = []
        for index in range(height):
            addition_polynomial = witt_sum_mod_p(characteristic, index)
            variables = addition_polynomial.parent().gens()
            x_variables = variables[: index + 1]
            y_variables = variables[index + 1 :]
            substitutions = {
                x_variables[degree]: z[degree] ** characteristic
                for degree in range(index + 1)
            }
            substitutions.update(
                {y_variables[degree]: -z[degree] for degree in range(index + 1)}
            )
            rhs.append(addition_polynomial.substitute(substitutions))
        return [-rhs[index] + z[index] ** characteristic - z[index] + f[index] for index in range(height)]

    def actions(z, _f, x, y):
        action = []
        for index in range(height):
            addition_polynomial = witt_sum_mod_p(characteristic, index)
            variables = addition_polynomial.parent().gens()
            x_variables = variables[: index + 1]
            y_variables = variables[index + 1 :]
            substitutions = {x_variables[degree]: z[degree] for degree in range(index + 1)}
            substitutions.update({y_variables[degree]: degree == 0 for degree in range(index + 1)})
            action.append(addition_polynomial.substitute(substitutions))
        return [action + [x, y]]

    return CoverTemplate(height, GF(characteristic), cover_group, equations, actions)


def witt_cover(functions, prec=10):
    '''Return the Artin-Schreier-Witt cover defined by quotient functions.'''
    from ascovers.as_covers.cover import ArtinSchreierCover

    quotient = functions[0].curve
    return ArtinSchreierCover(
        quotient,
        witt_template(quotient.characteristic, len(functions)),
        functions,
        branch_points=[],
        prec=prec,
    )


def quaternion_template():
    '''Return the characteristic-two quaternion template.'''
    cover_group = quaternion_group()

    def equations(z, f, _x, _y):
        return [f[0], f[1], f[2] + z[0] * f[0] + z[1] * (f[0] + f[1])]

    def actions(z, _f, x, y):
        return [
            [z[0] + 1, z[1], z[2] + z[0], x, y],
            [z[0], z[1] + 1, z[2] + z[1] + z[0], x, y],
        ]

    return CoverTemplate(3, GF(2), cover_group, equations, actions)


def quaternion_cover(functions, prec=10):
    '''Return the quaternion cover defined by three quotient functions.'''
    from ascovers.as_covers.cover import ArtinSchreierCover

    return ArtinSchreierCover(
        functions[0].curve,
        quaternion_template(),
        functions,
        branch_points=[],
        prec=prec,
    )


def hypoelementary_template(characteristic, order, character_value, zeta):
    '''Return the one-generator hypoelementary template.

    This reproduces the unfinished legacy template, but fixes the shape of the
    generator-action list so each generator has one action row.
    '''
    cover_group = hypoelementary_group(characteristic, order, character_value)

    def equations(z, f, _x, _y):
        return [
            1 / (zeta - character_value) * f[0] ** characteristic * z[0] ** characteristic
            - 1 / (zeta - character_value) * f[0] * z[0]
        ]

    def actions(z, f, x, y):
        return [
            [character_value * z[0] + f[0] * y, x, zeta * y],
            [z[0] + 1, x, y],
        ]

    return CoverTemplate(1, GF(characteristic), cover_group, equations, actions)


def dihedral_template(characteristic, height):
    '''Return the dihedral template based on Witt coordinates.'''
    base_template = witt_template(characteristic, height)
    cover_group = dihedral_group(characteristic, height)
    _field, z, _f, x, y = base_template.fct_field
    symmetry = [-z[index] for index in range(height)] + [x, -y]
    return CoverTemplate(
        height,
        GF(characteristic),
        cover_group,
        list(base_template.fcts),
        list(base_template.gp_action) + [symmetry],
    )


def dihedral_cover(functions, prec=10):
    '''Return the dihedral cover defined by quotient functions.'''
    from ascovers.as_covers.cover import ArtinSchreierCover

    quotient = functions[0].curve
    return ArtinSchreierCover(
        quotient,
        dihedral_template(quotient.characteristic, len(functions)),
        functions,
        branch_points=[],
        prec=prec,
    )


def d8_template():
    '''Return the characteristic-two D8 template from the legacy code.'''
    cover_group = dihedral_group(2, 2)

    def equations(z, f, _x, _y):
        return [f[0], f[1], f[1] * z[0]]

    def actions(z, _f, x, y):
        return [
            [z[0] + 1, z[1] + 1, z[2] + z[1] + 1, x, y],
            [z[0], z[1] + 1, z[2], x, y],
        ]

    return CoverTemplate(3, GF(2), cover_group, equations, actions)


def d8_cover(functions, prec=10):
    '''Return the D8 cover defined by quotient functions.'''
    from ascovers.as_covers.cover import ArtinSchreierCover

    return ArtinSchreierCover(
        functions[0].curve,
        d8_template(),
        functions,
        branch_points=[],
        prec=prec,
    )

def A4_template():
    '''Return the A4 template.'''
    cover_group = A4_group()
    F = GF(4)
    zeta = F.gens()[0]
    def equations(z, f, _x, _y):
        return [f[0], zeta*f[0]]

    def actions(z, _f, x, y):
        return [
            [z[0], z[1] + 1, x, y],
            [z[1], z[1] + z[0], zeta*x, zeta*y],
        ]

    return CoverTemplate(2, F, cover_group, equations, actions)


def A4_cover(function, prec=10):
    '''Return the A4 cover defined by quotient functions.'''
    from ascovers.as_covers.cover import ArtinSchreierCover
    polynomial = function.function.numerator()
    F = polynomial.parent().base_ring()
    x, y = polynomial.parent().gens()
    Rxy = PolynomialRing(F, 'x, y')
    Rx = PolynomialRing(F, 'x')
    polynomial = Rx(polynomial(x=x, y = 1))
    zeta = F.gens()[0] #this fails, if F != GF(4)

    if sum([coeff for index, coeff in enumerate(list(polynomial)) if index%3 == 0 ]) != 0:
        raise ValueError("The polynomial needs to satisfy tr(f) = 0")

    quotient = function.curve

    if quotient.exponent != 1 or quotient.polynomial.degree() != 1:
        raise ValueError("Only covers of P1 are implemented")

    return ArtinSchreierCover(
        quotient,
        A4_template(),
        [function, zeta * function],
        branch_points=[],
        prec=prec,
    )

template = CoverTemplate
witt_pol = witt_polynomial
dihedreal_template = dihedral_template
dihedreal_cover = dihedral_cover
D8_template = d8_template
D8_cover = d8_cover
