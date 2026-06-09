# ascovers

`ascovers` is a SageMath package for superelliptic curves and their
`p`-group covers in characteristic `p`.

The migrated package lives under `src/ascovers`.  The legacy Sage scripts are
kept in `old/` while the package is rewritten module by module.

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
`superelliptic_form`.
