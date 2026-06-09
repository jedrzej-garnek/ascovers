'''Superelliptic curve functionality.'''

from ascovers.superelliptic.curve import (
    SuperellipticCurve,
    reduction,
    reduction_form,
    superelliptic,
)
from ascovers.superelliptic.exceptions import PendingMigrationError
from ascovers.superelliptic.form import SuperellipticForm, superelliptic_form
from ascovers.superelliptic.function import SuperellipticFunction, superelliptic_function

__all__ = [
    "PendingMigrationError",
    "SuperellipticCurve",
    "SuperellipticForm",
    "SuperellipticFunction",
    "reduction",
    "reduction_form",
    "superelliptic",
    "superelliptic_form",
    "superelliptic_function",
]
