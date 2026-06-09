# ascovers

`ascovers` is a SageMath package for superelliptic curves and their
`p`-group covers in characteristic `p`.

The migrated package lives under `src/ascovers`.  The legacy Sage scripts are
kept in `old/` while the package is rewritten module by module.

## Documentation

- [Superelliptic curves](docs/superelliptic.md)

## Quick Start

Use the package from a SageMath Python environment:

```python
from sage.all import GF, PolynomialRing
from ascovers import SuperellipticCurve

F = GF(7)
R = PolynomialRing(F, "x")
x = R.gen()

C = SuperellipticCurve(x**5 + x, 4)
print(C.genus())
print(C.holomorphic_differentials_basis())
```

## Development

Run the tests from a SageMath Python environment:

```bash
conda activate sage
PYTHONPATH=src sage -python -m unittest discover -s test -v
```

If Sage cannot write to its default dot-directory, set `DOT_SAGE` to a
writable directory, for example `DOT_SAGE=/tmp/ascovers-sage-dot`.

The first migrated module is `ascovers.superelliptic.curve`, which defines
`SuperellipticCurve` and the compatibility alias `superelliptic`.  The
function and form modules define `SuperellipticFunction`,
`SuperellipticForm`, and the legacy aliases `superelliptic_function` and
`superelliptic_form`.  Cech-de Rham cocycles are represented by
`SuperellipticDeRhamCocycle`, with legacy alias `superelliptic_cech`.
