'''Compatibility exports for auxiliary AS-cover computations.'''

from ascovers.as_covers.matrices import (
    magma_is_isomorphic,
    magma_module_decomposition,
    magma_module_decomposition3,
)

__all__ = [
    "magma_is_isomorphic",
    "magma_module_decomposition",
    "magma_module_decomposition3",
]
