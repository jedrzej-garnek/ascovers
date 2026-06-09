'''Artin-Schreier type p-group covers of superelliptic curves.'''

from itertools import product
from random import randrange

from sage.all import Integer, LaurentSeriesRing, PolynomialRing, gcd, matrix, prod, vector, xgcd

from ascovers.as_covers.holomorphic import (
    holomorphic_combinations,
    holomorphic_combinations_forms,
    holomorphic_combinations_functions,
    holomorphic_combinations_mixed,
)
from ascovers.as_covers.matrices import group_action_matrices
from ascovers.as_covers.transform import artin_schreier_transform


class ArtinSchreierCover:
    '''A p-group cover of a superelliptic curve defined by an AS template.'''

    def __init__(self, quotient, cover_template, functions, branch_points=None, prec=10):
        '''Create the cover obtained by substituting quotient functions into a template.

        INPUT:

        - ``quotient`` -- a ``SuperellipticCurve``.
        - ``cover_template`` -- a ``CoverTemplate`` describing equations/actions.
        - ``functions`` -- quotient functions substituted for the template ``f_i``.
        - ``branch_points`` -- additional branch points; points at infinity are included.
        - ``prec`` -- Laurent precision used for local expansions.
        '''
        from ascovers.as_covers.form import ArtinSchreierForm
        from ascovers.as_covers.function import ArtinSchreierFunction

        if branch_points is None:
            branch_points = []

        self.quotient = quotient
        self.cover_template = cover_template
        self.functions = [
            function if getattr(function, "curve", None) is quotient else quotient.function(function)
            for function in functions
        ]
        self.height = int(cover_template.height)
        if len(self.functions) != self.height:
            raise ValueError("the number of quotient functions must match the template height")

        self.base_ring = quotient.base_ring
        self.characteristic = quotient.characteristic
        self.prec = int(prec)
        self.group = cover_template.group
        self.nb_of_pts_at_infty = gcd(quotient.exponent, quotient.polynomial.degree())
        self.branch_points = list(range(int(self.nb_of_pts_at_infty))) + list(branch_points)

        self.x_series = {}
        self.y_series = {}
        self.z_series = {}
        self.dx_series = {}
        self.jumps = {}
        self._compute_local_expansions()

        variable_names = ["x", "y"] + [f"z{index}" for index in range(self.height)]
        polynomial_ring = PolynomialRing(self.base_ring, tuple(variable_names))
        function_field = polynomial_ring.fraction_field()
        x_coordinate, y_coordinate = polynomial_ring.gens()[:2]
        z_coordinates = polynomial_ring.gens()[2:]
        self.fct_field = (
            function_field,
            polynomial_ring,
            x_coordinate,
            y_coordinate,
            z_coordinates,
        )

        self.x = ArtinSchreierFunction(self, x_coordinate)
        self.y = ArtinSchreierFunction(self, y_coordinate)
        self.z = [ArtinSchreierFunction(self, z_coordinate) for z_coordinate in z_coordinates]
        self.dx = ArtinSchreierForm(self, 1)
        self.one = ArtinSchreierFunction(self, 1)
        self.rhs = self._substitute_template_equations()
        self.dz = self._differentiate_cover_equations()

    def __repr__(self):
        '''Return a readable description of the cover equations.'''
        result = f"({self.group.short_name})-cover of {self.quotient} with equations:\n"
        for index in range(self.height):
            rhs = str(self.cover_template.fcts[index])
            for function_index, function in enumerate(self.functions):
                rhs = rhs.replace(f"f{function_index}", f"({function})")
            result += f"z{index}^p - z{index} = {rhs}\n"
        return result

    def _quotient_expansion(self, function, point):
        '''Return the expansion of a quotient function at a supported point.'''
        if isinstance(point, (int, Integer)):
            return function.expansion_at_infty(place=int(point), prec=self.prec)
        return function.expansion(point, prec=self.prec)

    def _compute_local_expansions(self):
        '''Compute expansions of x, y, z_i, and dx at each tracked branch point.'''
        template_field, template_z, template_f, template_x, template_y = self.cover_template.fct_field

        for point in self.branch_points:
            x_series = self._quotient_expansion(self.quotient.x, point)
            y_series = self._quotient_expansion(self.quotient.y, point)
            quotient_function_series = [
                self._quotient_expansion(function, point) for function in self.functions
            ]
            z_series = []
            jumps = []

            for index in range(self.height):
                substitutions = {
                    template_z[j]: z_series[j] for j in range(index)
                }
                substitutions.update(
                    {template_z[j]: 0 for j in range(index, self.height)}
                )
                substitutions.update(
                    {
                        template_f[j]: quotient_function_series[j]
                        for j in range(self.height)
                    }
                )
                substitutions.update({template_x: x_series, template_y: y_series})
                local_rhs = template_field(self.cover_template.fcts[index]).substitute(substitutions)
                jump, _correction, old_parameter, generator_series = artin_schreier_transform(
                    local_rhs,
                    prec=self.prec,
                )
                x_series = x_series(old_parameter)
                y_series = y_series(old_parameter)
                z_series = [series(old_parameter) for series in z_series]
                quotient_function_series = [
                    series(old_parameter) for series in quotient_function_series
                ]
                z_series.append(generator_series)
                jumps.append(jump)

            self.jumps[point] = jumps
            self.x_series[point] = x_series
            self.y_series[point] = y_series
            self.z_series[point] = z_series
            self.dx_series[point] = x_series.derivative()

    def _lift_quotient_expression(self, expression):
        '''Coerce an expression in quotient coordinates into the cover function field.'''
        cover_field, cover_ring, x_coordinate, y_coordinate, _z_coordinates = self.fct_field
        quotient_field, _quotient_ring, quotient_x, quotient_y = self.quotient.fct_field
        quotient_expression = quotient_field(expression)
        numerator = quotient_expression.numerator().substitute(
            {quotient_x: x_coordinate, quotient_y: y_coordinate}
        )
        denominator = quotient_expression.denominator().substitute(
            {quotient_x: x_coordinate, quotient_y: y_coordinate}
        )
        return cover_field(cover_ring(numerator)) / cover_field(cover_ring(denominator))

    def _substitute_template_equations(self):
        '''Substitute concrete quotient functions into the template equations.'''
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        template_field, template_z, template_f, template_x, template_y = self.cover_template.fct_field
        substitutions = {
            template_z[index]: function_field(z_coordinates[index])
            for index in range(self.height)
        }
        substitutions.update(
            {
                template_f[index]: self._lift_quotient_expression(self.functions[index].function)
                for index in range(self.height)
            }
        )
        substitutions.update({template_x: x_coordinate, template_y: y_coordinate})
        return [
            function_field(template_field(self.cover_template.fcts[index]).substitute(substitutions))
            for index in range(self.height)
        ]

    def _differentiate_cover_equations(self):
        '''Compute coefficients of ``dz_i`` from the differentiated equations.'''
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        dy_coefficient = self._lift_quotient_expression(self.quotient.y.diffn().form.function)
        differentials = []
        for index in range(self.height):
            rhs = self.rhs[index]
            derivative = rhs.derivative(x_coordinate)
            derivative += rhs.derivative(y_coordinate) * dy_coefficient
            for previous_index in range(index):
                derivative += rhs.derivative(z_coordinates[previous_index]) * differentials[previous_index]
            differentials.append(-function_field(derivative))
        return differentials

    def function(self, expression):
        '''Coerce an expression to a function on this cover.'''
        from ascovers.as_covers.function import ArtinSchreierFunction

        return ArtinSchreierFunction(self, expression)

    def form(self, expression):
        '''Coerce an expression to a differential form on this cover.'''
        from ascovers.as_covers.form import ArtinSchreierForm

        return ArtinSchreierForm(self, expression)

    def cech(self, omega, function):
        '''Build a Cech-de Rham representative on this cover.'''
        from ascovers.as_covers.cech import ArtinSchreierDeRhamCocycle

        return ArtinSchreierDeRhamCocycle(self, omega, function)

    def exponent_of_different(self, place=0):
        '''Return the exponent of the different at the selected branch point.'''
        values = [0]
        for index in range(1, self.height + 1):
            jump = self.jumps[place][index - 1]
            if jump == 0:
                values.append(values[index - 1])
            else:
                values.append((jump + 1) * (self.characteristic - 1) + self.characteristic * values[index - 1])
        return values[self.height]

    def exponent_of_different_prim(self, place=0):
        '''Return the variant of the different exponent using jumps rather than jumps+1.'''
        values = [0]
        for index in range(1, self.height + 1):
            jump = self.jumps[place][index - 1]
            if jump == 0:
                values.append(values[index - 1])
            else:
                values.append(jump * (self.characteristic - 1) + self.characteristic * values[index - 1])
        return values[self.height]

    def exponent_of_different_bis(self, place=0):
        '''Return the variant of the different exponent using jumps-1.'''
        values = [0]
        for index in range(1, self.height + 1):
            jump = self.jumps[place][index - 1]
            if jump == 0:
                values.append(values[index - 1])
            else:
                values.append((jump - 1) * (self.characteristic - 1) + self.characteristic * values[index - 1])
        return values[self.height]

    def genus(self):
        '''Return the genus computed by Riemann-Hurwitz from local jumps.'''
        quotient_genus = self.quotient.genus()
        group_order = self.characteristic ** self.height
        different_contribution = sum(
            self.exponent_of_different(place) * len(self.fiber(place))
            for place in self.branch_points
        )
        return group_order * (quotient_genus - 1) + 1 + different_contribution / 2

    def stabilizer(self, place=0):
        '''Return the stabilizer of a selected point above a branch point.'''
        result = []
        for element in self.group.elts:
            fixed = True
            for index in range(self.height):
                if (
                    self.jumps[place][index] == 0
                    and self.z[index].group_action(element).diffn() - self.z[index].diffn() == 0 * self.dx
                    and self.z[index].group_action(element) != self.z[index]
                ):
                    fixed = False
            if fixed:
                result.append(element)
        return result

    def fiber(self, place=0):
        '''Return representatives of ``G/G_P`` for the point over ``place``.'''
        result = [self.group.one]
        stabilizer = self.stabilizer(place=place)
        for element in self.group.elts:
            if element == self.group.one:
                continue
            is_new_coset = True
            for representative in result:
                difference = (self.group.elt(element) * (-self.group.elt(representative))).as_tuple
                if difference in stabilizer:
                    is_new_coset = False
            if is_new_coset:
                result.append(element)
        return result

    def holomorphic_differentials_basis(self, threshold=8):
        '''Search for a basis of holomorphic differentials on the cover.'''
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        quotient_exponent = int(self.quotient.exponent)
        quotient_degree = int(self.quotient.polynomial.degree())
        candidate_items = []
        z_exponent_ranges = [range(int(self.characteristic)) for _index in range(self.height)]

        for x_power in range(0, int(threshold) * quotient_degree):
            for y_denominator_power in range(0, quotient_exponent):
                for z_exponents in product(*z_exponent_ranges):
                    coefficient = x_coordinate ** x_power / y_coordinate ** y_denominator_power
                    coefficient *= prod(
                        z_coordinates[index] ** int(z_exponents[index])
                        for index in range(self.height)
                    )
                    form = self.form(function_field(coefficient))
                    candidate_items.append((form, form.expansion_at_infty(place=self.branch_points[0])))

        forms = holomorphic_combinations(candidate_items)
        for place in self.branch_points:
            for group_element in self.fiber(place=place):
                if place != self.branch_points[0] or group_element != self.group.one:
                    forms = [
                        (form, form.group_action(group_element).expansion_at_infty(place=place))
                        for form in forms
                    ]
                    forms = holomorphic_combinations(forms)

        if len(forms) < self.genus():
            raise ValueError(f"found only {len(forms)} of {self.genus()} holomorphic forms; increase threshold")
        if len(forms) > self.genus():
            raise ValueError("found too many forms; increase precision or reduce the candidate space")
        return [form.reduce() for form in forms]

    def cartier_matrix(self, prec=50):
        '''Return the Cartier-Manin matrix on holomorphic differentials.'''
        genus = int(self.genus())
        result = matrix(self.base_ring, genus, genus)
        basis = self.holomorphic_differentials_basis()
        for index, omega in enumerate(basis):
            result[:, index] = vector(omega.cartier().coordinates(basis=basis))
        return result

    def at_most_poles(self, pole_order, threshold=8):
        '''Search for functions with pole order at infinity at most ``pole_order``.'''
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        quotient_exponent = int(self.quotient.exponent)
        quotient_degree = int(self.quotient.polynomial.degree())
        z_exponent_ranges = [range(int(self.characteristic)) for _index in range(self.height)]
        candidate_items = []

        for x_power in range(0, int(threshold) * quotient_degree):
            for y_power in range(0, quotient_exponent):
                for z_exponents in product(*z_exponent_ranges):
                    expression = x_coordinate ** x_power * y_coordinate ** y_power
                    expression *= prod(
                        z_coordinates[index] ** int(z_exponents[index])
                        for index in range(self.height)
                    )
                    function = self.function(function_field(expression))
                    candidate_items.append((function, function.expansion_at_infty()))

        functions = holomorphic_combinations_functions(candidate_items, pole_order)
        for place in self.branch_points[1:]:
            functions = [
                (function, function.expansion_at_infty(place=place))
                for function in functions
            ]
            functions = holomorphic_combinations_functions(functions, pole_order)
        return functions

    def at_most_poles_forms(self, pole_order, threshold=8):
        '''Search for forms with pole order at infinity at most ``pole_order``.'''
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        quotient_exponent = int(self.quotient.exponent)
        quotient_degree = int(self.quotient.polynomial.degree())
        z_exponent_ranges = [range(int(self.characteristic)) for _index in range(self.height)]
        candidate_items = []

        for x_power in range(0, int(threshold) * quotient_degree):
            for y_denominator_power in range(0, quotient_exponent):
                for z_exponents in product(*z_exponent_ranges):
                    coefficient = x_coordinate ** x_power / y_coordinate ** y_denominator_power
                    coefficient *= prod(
                        z_coordinates[index] ** int(z_exponents[index])
                        for index in range(self.height)
                    )
                    form = self.form(function_field(coefficient))
                    candidate_items.append((form, form.expansion_at_infty()))

        forms = holomorphic_combinations_forms(candidate_items, pole_order)
        for place in self.branch_points[1:]:
            forms = [(form, form.expansion_at_infty(place=place)) for form in forms]
            forms = holomorphic_combinations_forms(forms, pole_order)
        return forms

    def magical_element(self, threshold=8):
        '''Return elements with bounded poles and nonzero trace.'''
        return [
            function
            for function in self.at_most_poles(self.exponent_of_different_prim(), threshold)
            if function.trace().function != 0
        ]

    def pseudo_magical_element(self, threshold=8):
        '''Return elements with different-bounded poles and nonzero trace.'''
        return [
            function
            for function in self.at_most_poles(self.exponent_of_different(), threshold)
            if function.trace().function != 0
        ]

    def uniformizer(self, place=0, threshold1=10, threshold2=10):
        '''Return a uniformizer at a selected point above infinity.'''
        bounded_functions = self.at_most_poles(threshold1, threshold=threshold2)
        simple_functions = self.z + [self.x, self.y]
        for first in bounded_functions:
            for second in simple_functions:
                divisor, first_exponent, second_exponent = xgcd(
                    first.valuation(place),
                    second.valuation(place),
                )
                if divisor == 1:
                    return first ** first_exponent * second ** second_exponent
        raise ValueError("increase the pole-search thresholds")

    def quasiuniformizer(self, place=0, threshold=10, auxiliary1=0, auxiliary2=0, **legacy_keywords):
        '''Return a function whose valuation is coprime to ``p``.

        ``auxiliary2`` is accepted for compatibility with the legacy signature;
        the old implementation did not use it in the final algorithm.
        '''
        if "auxilliary1" in legacy_keywords:
            auxiliary1 = legacy_keywords.pop("auxilliary1")
        if "auxilliary2" in legacy_keywords:
            auxiliary2 = legacy_keywords.pop("auxilliary2")
        if "auxiliary2" in legacy_keywords:
            auxiliary2 = legacy_keywords.pop("auxiliary2")
        if legacy_keywords:
            unknown = ", ".join(sorted(legacy_keywords))
            raise TypeError(f"unknown keyword argument(s): {unknown}")

        series_ring = LaurentSeriesRing(self.base_ring, "t", default_prec=self.prec)
        characteristic = self.characteristic
        if self.height == 1:
            return self.z[0]

        def difference_of_powers(first, second):
            first_series = series_ring(first.expansion_at_infty(place=place))
            second_series = series_ring(second.expansion_at_infty(place=place))
            first_valuation = first_series.valuation()
            second_valuation = second_series.valuation()
            common_divisor = gcd(first_valuation, second_valuation)
            if common_divisor:
                first_reduced = first_valuation // common_divisor
                second_reduced = second_valuation // common_divisor
            else:
                first_reduced, second_reduced = 1, 1
            if first_reduced < 0:
                first_reduced, second_reduced = -first_reduced, -second_reduced
            leading_ratio = list(first_series)[0] / list(second_series)[0]
            return first ** second_reduced - (leading_ratio * second) ** second_reduced

        def quotient_of_powers(first, second):
            first_series = series_ring(first.expansion_at_infty(place=place))
            second_series = series_ring(second.expansion_at_infty(place=place))
            first_valuation = first_series.valuation()
            second_valuation = second_series.valuation()
            common_divisor = gcd(first_valuation, second_valuation)
            if common_divisor:
                first_reduced = first_valuation // common_divisor
                second_reduced = second_valuation // common_divisor
            else:
                first_reduced, second_reduced = 1, 1
            leading_coefficient = list(first_series ** second_reduced / second_series ** first_reduced)[0]
            return first ** second_reduced / second ** first_reduced - leading_coefficient * self.one

        if auxiliary1 == 0:
            auxiliary1 = self.z[-1]

        candidate = quotient_of_powers(self.z[-2], auxiliary1)
        while candidate.valuation(place=place) % characteristic == 0:
            if randrange(2) == 0:
                candidate = difference_of_powers(auxiliary1, candidate)
            else:
                candidate = quotient_of_powers(auxiliary1, candidate)
        return candidate

    def ith_ramification_group(self, index, place=0, quasiuniformizer=0, threshold=20):
        '''Return the ``index``-th lower ramification group at a point above infinity.'''
        if isinstance(quasiuniformizer, (int, Integer)):
            uniformizer = self.quasiuniformizer(place, threshold=threshold)
        else:
            uniformizer = quasiuniformizer

        result = [self.group.elts[0]]
        for element in self.group.elts[1:]:
            translated = uniformizer.group_action(element)
            valuation = (translated - uniformizer).valuation(place)
            if valuation >= index + uniformizer.valuation(place):
                result.append(element)
        return result

    def ramification_jumps(self, place=0, quasiuniformizer=0, threshold=20):
        '''Return lower ramification jumps at a point above infinity.'''
        stabilizer = self.stabilizer(place=place)
        if isinstance(quasiuniformizer, (int, Integer)):
            uniformizer = self.quasiuniformizer(place, threshold=threshold)
        else:
            uniformizer = quasiuniformizer

        result = []
        index = 0
        current_group = list(stabilizer)
        while len(current_group) > 1:
            next_group = [current_group[0]]
            for element in current_group[1:]:
                translated = uniformizer.group_action(element)
                valuation = (translated - uniformizer).valuation(place)
                if valuation >= index + uniformizer.valuation(place):
                    next_group.append(element)
            if len(next_group) < len(current_group):
                result.append(index - 1)
            current_group = next_group
            index += 1
        return result

    def upper_ramification_jumps(self, place=0, quasiuniformizer=0, threshold=20):
        '''Return upper ramification jumps at a point above infinity.'''
        lower_jumps = self.ramification_jumps(
            place=place,
            quasiuniformizer=quasiuniformizer,
            threshold=threshold,
        )
        if not lower_jumps:
            return []
        result = [lower_jumps[0]]
        for index in range(1, len(lower_jumps)):
            denominator = len(self.stabilizer(place=place)) // len(
                self.ith_ramification_group(
                    lower_jumps[index],
                    place=place,
                    quasiuniformizer=quasiuniformizer,
                    threshold=threshold,
                )
            )
            result.append(result[index - 1] + (lower_jumps[index] - lower_jumps[index - 1]) // denominator)
        return result

    ith_ramification_gp = ith_ramification_group

    def a_number(self):
        '''Return the a-number computed from the Cartier matrix.'''
        return self.genus() - self.cartier_matrix().rank()

    def cohomology_of_structure_sheaf_basis(self, holo_basis=0, threshold=8):
        '''Return a basis of H^1(X, O_X) by Serre duality.'''
        if holo_basis == 0:
            holo_basis = self.holomorphic_differentials_basis(threshold=threshold)
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        quotient_exponent = int(self.quotient.exponent)
        z_exponent_ranges = [range(int(self.characteristic)) for _index in range(self.height)]
        vector_space = self.base_ring ** int(self.genus())
        span = vector_space.subspace([])
        result = []
        x_pole_power = 0
        while len(result) < self.genus():
            for y_denominator_power in range(0, quotient_exponent):
                for z_exponents in product(*z_exponent_ranges):
                    expression = prod(
                        z_coordinates[index] ** int(z_exponents[index])
                        for index in range(self.height)
                    ) / x_coordinate ** x_pole_power / y_coordinate ** y_denominator_power
                    function = self.function(function_field(expression))
                    pairings = [omega.serre_duality_pairing(function) for omega in holo_basis]
                    vector_pairing = vector(self.base_ring, pairings)
                    if vector_pairing not in span:
                        span = span + vector_space.subspace([vector_pairing])
                        result.append(function)
            x_pole_power += 1
        return result

    def lift_to_de_rham(self, function, basis=0, threshold=30):
        '''Find a form eta such that eta - d(function) is regular at infinity.'''
        if basis == 0:
            holomorphic_basis = self.holomorphic_differentials_basis(threshold=threshold)
        else:
            holomorphic_basis = basis

        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        quotient_exponent = int(self.quotient.exponent)
        quotient_degree = int(self.quotient.polynomial.degree())
        z_exponent_ranges = [range(int(self.characteristic)) for _index in range(self.height)]
        candidates = [(function.diffn(), function.diffn().expansion_at_infty())]
        for x_power in range(0, int(threshold) * quotient_degree):
            for y_denominator_power in range(0, quotient_exponent):
                for z_exponents in product(*z_exponent_ranges):
                    coefficient = x_coordinate ** x_power / y_coordinate ** y_denominator_power
                    coefficient *= prod(
                        z_coordinates[index] ** int(z_exponents[index])
                        for index in range(self.height)
                    )
                    form = self.form(function_field(coefficient))
                    candidates.append((form, form.expansion_at_infty()))
        forms = holomorphic_combinations(candidates)

        for place in self.branch_points:
            for group_element in self.fiber(place=place):
                if place != self.branch_points[0] or group_element != self.group.one:
                    forms = [
                        (form, form.group_action(group_element).expansion_at_infty(place=place))
                        for form in forms
                    ]
                    forms = holomorphic_combinations(forms)

        if len(forms) <= self.genus():
            raise ValueError("increase threshold")
        for form in forms:
            for scalar in self.base_ring:
                candidate = scalar * form + function.diffn()
                if candidate.is_regular_on_U0():
                    return candidate
        raise ValueError("could not find a de Rham lift")

    def lift_to_de_rham_form(self, form, threshold=20):
        '''Find Cech representatives cancelling a pole of the given form.'''
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        quotient_exponent = int(self.quotient.exponent)
        quotient_degree = int(self.quotient.polynomial.degree())
        z_exponent_ranges = [range(int(self.characteristic)) for _index in range(self.height)]
        candidates = [(form, form.expansion_at_infty())]

        for x_power in range(0, int(threshold) * quotient_degree):
            for y_denominator_power in range(0, quotient_exponent):
                for z_exponents in product(*z_exponent_ranges):
                    z_part = prod(
                        z_coordinates[index] ** int(z_exponents[index])
                        for index in range(self.height)
                    )
                    for expression in (
                        z_part / x_coordinate ** x_power / y_coordinate ** y_denominator_power,
                        z_part * x_coordinate ** x_power * y_coordinate ** y_denominator_power,
                    ):
                        function = self.function(function_field(expression))
                        exact_expansion = function.diffn().expansion_at_infty()
                        if exact_expansion != 0:
                            candidates.append((function, exact_expansion))

        combinations = holomorphic_combinations_mixed(candidates)
        if len(combinations) <= self.genus():
            raise ValueError("increase threshold")
        return [
            self.cech(form_part, -function_part)
            for form_part, function_part in combinations
            if form_part != 0 * self.dx
        ]

    def de_rham_basis(self, holo_basis=0, cohomology_basis=0, threshold=30):
        '''Return a Cech-de Rham basis built from holomorphic and Cech classes.'''
        if holo_basis == 0:
            holo_basis = self.holomorphic_differentials_basis(threshold=threshold)
        if cohomology_basis == 0:
            cohomology_basis = self.cohomology_of_structure_sheaf_basis(
                holo_basis=holo_basis,
                threshold=threshold,
            )
        result = [self.cech(omega, 0 * self.x) for omega in holo_basis]
        for function in cohomology_basis:
            result.append(self.cech(self.lift_to_de_rham(function, basis=holo_basis, threshold=threshold), function))
        return result

    def holo_polydifferentials_basis(self, mult, threshold=8):
        '''Return a basis of holomorphic polydifferentials of tensor degree ``mult``.'''
        from ascovers.as_covers.polyforms import holomorphic_polydifferentials_basis

        return holomorphic_polydifferentials_basis(self, mult, threshold=threshold)

    holomorphic_polydifferentials_basis = holo_polydifferentials_basis

    def symmetric_power_basis(self, tensor_power, threshold=8):
        '''Return the natural basis of a symmetric power of holomorphic forms.'''
        from ascovers.as_covers.polyforms import symmetric_power_basis

        return symmetric_power_basis(self, tensor_power, threshold=threshold)

    def canonical_ideal(self, tensor_power, threshold=8):
        '''Return the degree-``tensor_power`` component of the canonical ideal.'''
        from ascovers.as_covers.polyforms import canonical_ideal

        return canonical_ideal(self, tensor_power, threshold=threshold)

    def canonical_ideal_polynomials(self, tensor_power, threshold=8):
        '''Return polynomial equations for a canonical ideal component.'''
        from ascovers.as_covers.polyforms import canonical_ideal_polynomials

        return canonical_ideal_polynomials(self, tensor_power, threshold=threshold)

    def group_action_canonical_ideal(self, tensor_power, threshold=8):
        '''Return group-action matrices on a canonical ideal component.'''
        from ascovers.as_covers.polyforms import group_action_canonical_ideal

        return group_action_canonical_ideal(self, tensor_power, threshold=threshold)

    def group_action_matrices_holo(self, basis=0, threshold=10):
        '''Return group-action matrices on holomorphic differentials.'''
        if basis == 0:
            basis = self.holomorphic_differentials_basis(threshold=threshold)
        return group_action_matrices(self.base_ring, basis, self.group.gens, basis=basis)

    def group_action_matrices_poly(self, mult, basis=0, threshold=10):
        '''Return group-action matrices on polydifferentials.'''
        if basis == 0:
            basis = self.holo_polydifferentials_basis(mult=mult, threshold=threshold)
        return group_action_matrices(self.base_ring, basis, self.group.gens, basis=basis)

    def group_action_matrices_dR(self, basis=0, threshold=8):
        '''Return group-action matrices on de Rham cohomology.'''
        if basis == 0:
            holomorphic_basis = self.holomorphic_differentials_basis(threshold=threshold)
            structure_basis = self.cohomology_of_structure_sheaf_basis(
                holo_basis=holomorphic_basis,
                threshold=threshold,
            )
            de_rham_basis = self.de_rham_basis(
                holo_basis=holomorphic_basis,
                cohomology_basis=structure_basis,
                threshold=threshold,
            )
            basis = [holomorphic_basis, structure_basis, de_rham_basis]
        return group_action_matrices(self.base_ring, basis[2], self.group.gens, basis=basis)

    def dR_pairing_matrix(self, basis=0, threshold=8):
        '''Return the de Rham pairing matrix in the selected basis.'''
        if basis == 0:
            holomorphic_basis = self.holomorphic_differentials_basis(threshold=threshold)
            structure_basis = self.cohomology_of_structure_sheaf_basis(
                holo_basis=holomorphic_basis,
                threshold=threshold,
            )
            de_rham_basis = self.de_rham_basis(
                holo_basis=holomorphic_basis,
                cohomology_basis=structure_basis,
                threshold=threshold,
            )
        else:
            holomorphic_basis, _structure_basis, de_rham_basis = basis
        genus = len(holomorphic_basis)
        result = matrix(self.base_ring, 2 * genus, 2 * genus)
        for row in range(2 * genus):
            for column in range(2 * genus):
                result[row, column] = de_rham_basis[row].pairing(de_rham_basis[column])
        return result

    def group_action_matrices_log_holo(self):
        '''Return group-action matrices on forms with at most logarithmic poles.'''
        basis = self.at_most_poles_forms(1)
        return group_action_matrices(self.base_ring, basis, self.group.gens, basis=basis)

    def riemann_roch_space(self, pole_orders, threshold=8):
        '''Find functions with pole bounds indexed by ``(place, group_element)``.'''
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        quotient_exponent = int(self.quotient.exponent)
        quotient_degree = int(self.quotient.polynomial.degree())
        z_exponent_ranges = [range(int(self.characteristic)) for _index in range(self.height)]
        candidate_items = []
        for x_power in range(0, int(threshold) * quotient_degree):
            for y_power in range(0, quotient_exponent):
                for z_exponents in product(*z_exponent_ranges):
                    expression = x_coordinate ** x_power * y_coordinate ** y_power
                    expression *= prod(
                        z_coordinates[index] ** int(z_exponents[index])
                        for index in range(self.height)
                    )
                    function = self.function(function_field(expression))
                    candidate_items.append((function, function.expansion_at_infty()))

        functions = holomorphic_combinations_functions(
            candidate_items,
            pole_orders[(0, self.group.one)],
        )
        for place in range(int(self.nb_of_pts_at_infty)):
            for group_element in self.fiber(place=place):
                if place != 0 or group_element != self.group.one:
                    functions = [
                        (
                            function,
                            function.group_action(group_element).expansion_at_infty(place=place),
                        )
                        for function in functions
                    ]
                    functions = holomorphic_combinations_functions(
                        functions,
                        pole_orders[(place, group_element)],
                    )
        return functions

    def riemann_roch_space_forms(self, pole_orders, threshold=8):
        '''Find forms with pole bounds indexed by ``(place, group_element)``.'''
        function_field, _polynomial_ring, x_coordinate, y_coordinate, z_coordinates = self.fct_field
        quotient_exponent = int(self.quotient.exponent)
        quotient_degree = int(self.quotient.polynomial.degree())
        z_exponent_ranges = [range(int(self.characteristic)) for _index in range(self.height)]
        candidate_items = []
        for x_power in range(0, int(threshold) * quotient_degree):
            for y_denominator_power in range(0, quotient_exponent):
                for z_exponents in product(*z_exponent_ranges):
                    coefficient = x_coordinate ** x_power / y_coordinate ** y_denominator_power
                    coefficient *= prod(
                        z_coordinates[index] ** int(z_exponents[index])
                        for index in range(self.height)
                    )
                    form = self.form(function_field(coefficient))
                    candidate_items.append((form, form.expansion_at_infty()))

        forms = holomorphic_combinations_forms(
            candidate_items,
            pole_orders[(0, self.group.one)],
        )
        for place in range(int(self.nb_of_pts_at_infty)):
            for group_element in self.fiber(place=place):
                if place != 0 or group_element != self.group.one:
                    forms = [
                        (
                            form,
                            form.group_action(group_element).expansion_at_infty(place=place),
                        )
                        for form in forms
                    ]
                    forms = holomorphic_combinations_forms(
                        forms,
                        pole_orders[(place, group_element)],
                    )
        return forms


as_cover = ArtinSchreierCover
