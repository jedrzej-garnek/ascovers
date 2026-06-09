'''Linear algebra helpers for imposing regularity from Laurent expansions.'''

from sage.all import FractionField, LaurentSeriesRing, matrix


def _series_ring(cover):
    '''Return the Laurent series ring used for expansions on an AS cover.'''
    return LaurentSeriesRing(cover.base_ring, "t", default_prec=cover.prec)


def _negative_coefficients(series, minimal_valuation, pole_order=0):
    '''Return coefficients that must vanish to impose the requested pole bound.'''
    series = series.parent()(series)
    if series == 0:
        return [0] * int(-minimal_valuation - pole_order)
    shift = int(-minimal_valuation + series.valuation())
    coefficients = [0] * shift + list(series) + [0] * int(-minimal_valuation)
    return coefficients[: int(-minimal_valuation - pole_order)]


def holomorphic_combinations(items):
    '''Return linear combinations of forms that are regular at a chosen place.

    INPUT:

    - ``items`` -- list of pairs ``(form, expansion)``.

    The expansions should all be taken at the same place.  The result is a basis
    of the kernel of the principal-part coefficient matrix.
    '''
    cover = items[0][0].curve
    base_field = cover.base_ring
    minimal_valuation = min(expansion.valuation() for _form, expansion in items)
    if minimal_valuation >= 0:
        return [form for form, _expansion in items]

    coefficient_rows = [
        _negative_coefficients(expansion, minimal_valuation, pole_order=0)
        for _form, expansion in items
    ]
    kernel = matrix(base_field, coefficient_rows).kernel()

    result = []
    for vector in kernel.basis():
        combination = 0 * items[0][0]
        for coefficient, (form, _expansion) in zip(vector, items):
            combination += coefficient * form
        result.append(combination)
    return result


def holomorphic_combinations_mixed(items):
    '''Return regular combinations of forms and exact differentials.

    Each item is ``(object, expansion)``, where the object is either an AS form or
    an AS function.  The returned pairs are ``(form_part, function_part)``.
    '''
    from ascovers.as_covers.form import ArtinSchreierForm
    from ascovers.as_covers.function import ArtinSchreierFunction

    cover = items[0][0].curve
    base_field = cover.base_ring
    minimal_valuation = min(expansion.valuation() for _object, expansion in items)
    if minimal_valuation >= 0:
        result = []
        for item, _expansion in items:
            if isinstance(item, ArtinSchreierForm):
                result.append((item, 0 * cover.x))
            elif isinstance(item, ArtinSchreierFunction):
                result.append((0 * cover.dx, item))
            else:
                raise TypeError("mixed regularity items must be AS forms or AS functions")
        return result

    coefficient_rows = [
        _negative_coefficients(expansion, minimal_valuation, pole_order=0)
        for _object, expansion in items
    ]
    kernel = matrix(base_field, coefficient_rows).kernel()

    result = []
    for vector in kernel.basis():
        form_part = 0 * cover.dx
        function_part = 0 * cover.x
        for coefficient, (item, _expansion) in zip(vector, items):
            if isinstance(item, ArtinSchreierForm):
                form_part += coefficient * item
            elif isinstance(item, ArtinSchreierFunction):
                function_part += coefficient * item
            else:
                raise TypeError("mixed regularity items must be AS forms or AS functions")
        result.append((form_part, function_part))
    return result


def holomorphic_combinations_functions(items, pole_order):
    '''Return combinations of functions with pole order at most ``pole_order``.'''
    cover = items[0][0].curve
    base_field = cover.base_ring
    series_ring = _series_ring(cover)
    minimal_valuation = min(series_ring(expansion).valuation() for _function, expansion in items)
    if minimal_valuation >= -pole_order:
        return [function for function, _expansion in items]

    coefficient_rows = [
        _negative_coefficients(series_ring(expansion), minimal_valuation, pole_order=pole_order)
        for _function, expansion in items
    ]
    kernel = matrix(base_field, coefficient_rows).kernel()

    result = []
    for vector in kernel.basis():
        combination = 0 * cover.x
        for coefficient, (function, _expansion) in zip(vector, items):
            combination += coefficient * function
        result.append(combination)
    return result


def holomorphic_combinations_forms(items, pole_order):
    '''Return combinations of forms with pole order at most ``pole_order``.'''
    cover = items[0][0].curve
    base_field = cover.base_ring
    series_ring = _series_ring(cover)
    minimal_valuation = min(series_ring(expansion).valuation() for _form, expansion in items)
    if minimal_valuation >= -pole_order:
        return [form for form, _expansion in items]

    coefficient_rows = [
        _negative_coefficients(series_ring(expansion), minimal_valuation, pole_order=pole_order)
        for _form, expansion in items
    ]
    kernel = matrix(base_field, coefficient_rows).kernel()

    result = []
    for vector in kernel.basis():
        combination = 0 * cover.dx
        for coefficient, (form, _expansion) in zip(vector, items):
            combination += coefficient * form
        result.append(combination)
    return result


holomorphic_combinations_fcts = holomorphic_combinations_functions
