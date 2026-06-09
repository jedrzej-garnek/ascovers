'''Matrix and Magma helpers for group actions on AS-cover spaces.'''

from sage.all import matrix, vector


def group_action_matrices(base_field, space, group_elements, basis, info=False):
    '''Return matrices for the action of selected group elements on a space.

    INPUT:

    - ``base_field`` -- coefficient field for the matrices.
    - ``space`` -- list of objects with ``group_action`` and ``coordinates``.
    - ``group_elements`` -- raw group elements whose actions are represented.
    - ``basis`` -- basis passed to ``coordinates``.
    - ``info`` -- print progress messages when true.
    '''
    dimension = len(space)
    matrices = [matrix(base_field, dimension, dimension) for _element in group_elements]
    for matrix_index, group_element in enumerate(group_elements):
        if info:
            print(f"Matrix for group element {group_element}")
        for column, element in enumerate(space):
            if info:
                print(f"coordinates of element {column + 1} out of {dimension}")
            transformed = element.group_action(group_element)
            matrices[matrix_index][:, column] = vector(base_field, transformed.coordinates(basis=basis))
    return matrices


def _matrix_entries_for_magma(action_matrix):
    '''Return matrix entries in Magma's flat matrix-constructor syntax.'''
    return str(list(action_matrix)).replace("(", "").replace(")", "")


def _magma_prefix(base_field, dimension):
    '''Return Magma declarations for a finite field and matrix algebra.'''
    field_order = base_field.order()
    characteristic = field_order.factor()[0][0]
    result = ""
    if field_order != characteristic:
        result += f"F<a> := GF({field_order});"
    return result + f"A := MatrixAlgebra<GF({field_order}),{dimension}|"


def magma_module_decomposition(first_matrix, second_matrix, text=False, prefix="", suffix="", matrices=True):
    '''Return or run Magma code decomposing a two-generator module.

    If ``text`` is true, the Magma command string is returned.  Otherwise this
    calls Sage's ``magma_free`` interface when available.
    '''
    base_field = first_matrix.parent().base_ring()
    dimension = first_matrix.dimensions()[0]
    command = prefix
    command += _magma_prefix(base_field, dimension)
    command += _matrix_entries_for_magma(first_matrix)
    command += ","
    command += _matrix_entries_for_magma(second_matrix)
    command += ">;"
    command += f"M := RModule(RSpace(GF({base_field.order()}),{dimension}), A);"
    command += "L := IndecomposableSummands(M); L;"
    if matrices:
        command += "for i in [1 .. #L] do print(Generators(Action(L[i]))); end for;"
    command += suffix
    if text:
        return command

    try:
        from sage.all import magma_free
    except ImportError as error:
        raise RuntimeError("Sage's magma_free interface is not available; call with text=True") from error
    return magma_free(command)


def magma_module_decomposition3(
    first_matrix,
    second_matrix,
    third_matrix,
    text=False,
    prefix="",
    suffix="",
    matrices=True,
):
    '''Return or run Magma code decomposing a three-generator module.'''
    base_field = first_matrix.parent().base_ring()
    dimension = first_matrix.dimensions()[0]
    command = prefix
    command += _magma_prefix(base_field, dimension)
    command += _matrix_entries_for_magma(first_matrix)
    command += ","
    command += _matrix_entries_for_magma(second_matrix)
    command += ","
    command += _matrix_entries_for_magma(third_matrix)
    command += ">;"
    command += f"M := RModule(RSpace(GF({base_field.order()}),{dimension}), A);"
    command += "L := IndecomposableSummands(M); L;"
    if matrices:
        command += "for i in [1 .. #L] do print(Generators(Action(L[i]))); end for;"
    command += suffix
    if text:
        return command

    try:
        from sage.all import magma_free
    except ImportError as error:
        raise RuntimeError("Sage's magma_free interface is not available; call with text=True") from error
    return magma_free(command)


def magma_is_isomorphic(first_a, first_b, second_a, second_b, text=False):
    '''Return or run Magma code testing isomorphism of two two-generator modules.'''
    base_field = first_a.parent().base_ring()
    dimension = first_a.dimensions()[0]
    command = ""
    field_order = base_field.order()
    characteristic = field_order.factor()[0][0]
    if field_order != characteristic:
        command += f"F<a> := GF({field_order});"
    command += f"A1 := MatrixAlgebra<GF({field_order}),{dimension}|"
    command += _matrix_entries_for_magma(first_a)
    command += ","
    command += _matrix_entries_for_magma(first_b)
    command += ">;"
    command += f"M1 := RModule(RSpace(GF({field_order}),{dimension}), A1);"
    command += f"A2 := MatrixAlgebra<GF({field_order}),{dimension}|"
    command += _matrix_entries_for_magma(second_a)
    command += ","
    command += _matrix_entries_for_magma(second_b)
    command += ">;"
    command += f"M2 := RModule(RSpace(GF({field_order}),{dimension}), A2);"
    command += "IsIsomorphic(M1, M2);"
    if text:
        return command

    try:
        from sage.all import magma_free
    except ImportError as error:
        raise RuntimeError("Sage's magma_free interface is not available; call with text=True") from error
    return magma_free(command)


as_group_action_matrices = group_action_matrices
