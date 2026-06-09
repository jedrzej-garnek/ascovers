'''Public package interface for ascovers.'''

from ascovers.superelliptic import (
    SuperellipticCurve,
    SuperellipticForm,
    SuperellipticFunction,
    reduction,
    reduction_form,
    superelliptic,
    superelliptic_form,
    superelliptic_function,
)

__all__ = [
    "SuperellipticCurve",
    "SuperellipticForm",
    "SuperellipticFunction",
    "reduction",
    "reduction_form",
    "superelliptic",
    "superelliptic_form",
    "superelliptic_function",
]

__version__ = "0.1.0"
