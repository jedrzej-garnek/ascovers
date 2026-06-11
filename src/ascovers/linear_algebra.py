'''Linear algebra helpers used by cohomology computations.'''

from sage.all import GF, Set, VectorSpace, block_matrix, matrix, vector


def linear_combination_coefficients(target_vector, basis_vectors):
    '''Return coefficients expressing a vector in a given spanning list.'''
    coefficient_matrix = block_matrix([[matrix([basis_vector])] for basis_vector in basis_vectors])
    coefficient_matrix = coefficient_matrix.transpose()
    return coefficient_matrix.solve_right(vector(target_vector))


def linear_combination_coeffs(target_vector, basis_vectors):
    '''Compatibility alias for the legacy helper name.'''
    return linear_combination_coefficients(target_vector, basis_vectors)


def from_coordinates(coordinates, basis):
    '''Construct a linear combination of basis elements from coordinates.'''
    result = 0 * basis[0]
    for index, coefficient in enumerate(coordinates):
        result += coefficient * basis[index]
    return result


def preimage(subspace, ambient_space, linear_map):
    '''Return the preimage of a subspace under a matrix.'''
    basis_preimage = linear_map.right_kernel().basis()
    image_intersection = subspace.intersection(linear_map.transpose().image())
    for vector_in_image in image_intersection.basis():
        basis_preimage.append(linear_map.solve_right(vector_in_image))
    return ambient_space.subspace(basis_preimage)


def image(subspace, ambient_space, linear_map):
    '''Return the image of a subspace under a matrix.'''
    basis_image = [linear_map * basis_vector for basis_vector in subspace.basis()]
    return ambient_space.subspace(basis_image)


def is_final(final_type, dimension):
    '''Return whether a sequence is a valid final type of the given dimension.'''
    if final_type[0] != 0:
        return False
    if final_type[-1] != dimension:
        return False
    for index in range(1, len(final_type)):
        previous = final_type[index - 1]
        current = final_type[index]
        if current != previous and current != previous + 1:
            return False
    return True


def flag(frobenius_matrix, verschiebung_matrix, characteristic, test=False):
    '''Compute the final type from Frobenius and Verschiebung matrices.'''
    dimension = frobenius_matrix.dimensions()[0]
    half_dimension = dimension // 2
    ambient_space = VectorSpace(GF(characteristic), dimension)

    flag_subspaces = [None] * (dimension + 1)
    active_dimensions = [False] * (dimension + 1)
    final_type = ["?"] * (dimension + 1)

    flag_subspaces[dimension] = ambient_space
    active_dimensions[dimension] = True

    while True in active_dimensions:
        index = active_dimensions.index(True)
        active_dimensions[index] = False
        subspace = flag_subspaces[index]

        image_subspace = image(subspace, ambient_space, verschiebung_matrix)
        image_dimension = image_subspace.dimension()
        final_type[index] = image_dimension

        preimage_subspace = preimage(subspace, ambient_space, frobenius_matrix)
        preimage_dimension = preimage_subspace.dimension()

        if flag_subspaces[image_dimension] is None:
            flag_subspaces[image_dimension] = image_subspace
            active_dimensions[image_dimension] = True

        if flag_subspaces[preimage_dimension] is None:
            flag_subspaces[preimage_dimension] = preimage_subspace
            active_dimensions[preimage_dimension] = True

    if test:
        print("(", final_type, ")")

    for index in range(dimension + 1):
        dual_index = dimension - index
        if final_type[index] == "?" and final_type[dual_index] != "?":
            final_type[index] = final_type[dual_index] - dual_index + half_dimension

    final_type[0] = 0

    for index in range(1, dimension + 1):
        if final_type[index] == "?":
            previous = final_type[index - 1]
            if previous != "?" and previous in final_type[index + 1:]:
                final_type[index] = previous

    for index in range(1, dimension + 1):
        if final_type[index] == "?":
            final_type[index] = min(final_type[index - 1] + 1, half_dimension)

    if is_final(final_type, half_dimension):
        return final_type[1:half_dimension + 1]

    raise ValueError(f"computed sequence is not a final type: {final_type[1:half_dimension + 1]}")


def jordan_block_dimensions(square_matrix):
    '''Return the dimensions of Jordan blocks of a square matrix.'''
    jordan_form = square_matrix.jordan_form()
    subdivisions = jordan_form.subdivisions()[0]
    matrix_dimension = square_matrix.dimensions()[0]
    if len(subdivisions) == 0:
        return [matrix_dimension]
    return (
        [subdivisions[0]]
        + [subdivisions[index + 1] - subdivisions[index] for index in range(len(subdivisions) - 1)]
        + [matrix_dimension - subdivisions[-1]]
    )


def jordan_dims(square_matrix):
    '''Compatibility alias for the legacy Jordan block helper.'''
    return jordan_block_dimensions(square_matrix)


def linear_representation_polynomials(polynomial, basis_polynomials):
    '''Return coefficients expressing a polynomial in a list of polynomials.'''
    base_ring = polynomial.parent().base_ring()
    monomials = list(polynomial.monomials())
    for basis_polynomial in basis_polynomials:
        monomials.extend(basis_polynomial.monomials())
    monomials = list(Set(monomials))

    coefficient_matrix = matrix(base_ring, len(monomials), len(basis_polynomials))
    for column_index, basis_polynomial in enumerate(basis_polynomials):
        for row_index, monomial in enumerate(monomials):
            coefficient_matrix[row_index, column_index] = basis_polynomial.monomial_coefficient(monomial)

    target_vector = vector([polynomial.monomial_coefficient(monomial) for monomial in monomials])
    return list(coefficient_matrix.solve_right(target_vector))
