'''Superelliptic curve functionality.'''

from ascovers.superelliptic.cech import (
    SuperellipticCech,
    SuperellipticCechCocycle,
    SuperellipticDeRhamCocycle,
    superelliptic_cech,
)
from ascovers.superelliptic.curve import (
    SuperellipticCurve,
    reduction,
    reduction_form,
    superelliptic,
)
from ascovers.superelliptic.decomposition import decomposition_g0_g8, decomposition_omega0_omega8
from ascovers.superelliptic.exceptions import PendingMigrationError
from ascovers.superelliptic.form import SuperellipticForm, superelliptic_form
from ascovers.superelliptic.function import SuperellipticFunction, superelliptic_function

__all__ = [
    "PendingMigrationError",
    "SuperellipticCech",
    "SuperellipticCechCocycle",
    "SuperellipticCurve",
    "SuperellipticDeRhamCocycle",
    "SuperellipticForm",
    "SuperellipticFunction",
    "reduction",
    "reduction_form",
    "decomposition_g0_g8",
    "decomposition_omega0_omega8",
    "superelliptic",
    "superelliptic_cech",
    "superelliptic_form",
    "superelliptic_function",
]
