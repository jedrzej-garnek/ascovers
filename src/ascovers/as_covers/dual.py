'''Trace-dual and normal-basis decomposition helpers for AS covers.'''

from sage.all import FractionField, matrix, vector


def trace_dual_element(cover, normal_basis_element):
    '''Return the trace-dual element to ``normal_basis_element``.

    The result ``z_dual`` satisfies
    ``Tr(normal_basis_element * g(z_dual)) = 1`` for the identity element and
    ``0`` for the other group elements, provided the normal-basis element is
    nondegenerate.
    '''
    from ascovers.as_covers.function import ArtinSchreierFunction

    function_field, polynomial_ring, _x_coordinate, _y_coordinate, _z_coordinates = cover.fct_field
    rational_field = FractionField(polynomial_ring)
    group_elements = list(cover.group.elts)
    dimension = len(group_elements)

    trace_matrix = matrix(rational_field, dimension, dimension)
    translates = [normal_basis_element.group_action(element) for element in group_elements]
    for row in range(dimension):
        for column in range(dimension):
            trace_value = (translates[row] * translates[column]).trace(super=False).function
            trace_matrix[row, column] = rational_field(trace_value)

    determinant = trace_matrix.determinant()
    if determinant == 0:
        raise ValueError("the element does not generate a nondegenerate trace pairing")

    result = ArtinSchreierFunction(cover, 0)
    identity_column = vector(
        rational_field,
        [rational_field(group_element == cover.group.one) for group_element in group_elements],
    )
    for column in range(dimension):
        auxiliary_matrix = matrix(rational_field, trace_matrix)
        auxiliary_matrix[:, column] = identity_column
        coefficient = auxiliary_matrix.determinant() / determinant
        result += ArtinSchreierFunction(cover, function_field(coefficient) * translates[column].function)
    return result.reduce()


def ith_magical_component(omega, dual_element, group_element, super=True):
    '''Return one component in a normal-basis decomposition of a form.'''
    translated_dual = dual_element.group_action(group_element)
    component = translated_dual * omega
    return component.trace(super=super)


def combination_components(omega, normal_basis_element, target_element):
    '''Recombine normal-basis components of ``omega`` with ``target_element``.

    If ``omega = sum_g g(normal_basis_element) omega_g``, this returns
    ``sum_g g(target_element) omega_g``.
    '''
    cover = omega.curve
    dual = trace_dual_element(cover, normal_basis_element)
    result = 0 * cover.dx
    for group_element in cover.group.elts:
        component = ith_magical_component(omega, dual, group_element, super=False)
        result += target_element.group_action(group_element) * component
    return result


dual_elt = trace_dual_element
