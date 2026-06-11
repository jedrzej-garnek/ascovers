# Artin-Schreier Type Covers

The `ascovers.as_covers` package handles `p`-group covers of superelliptic
curves in characteristic `p`.  A cover is built from:

- a quotient `SuperellipticCurve`,
- a finite group model,
- a symbolic `CoverTemplate` describing equations and generator actions,
- quotient functions substituted into the template parameters.

The main class is `ArtinSchreierCover`.  The legacy aliases
`as_cover`, `as_function`, `as_form`, `as_cech`, `group`, and `template` remain
available for older code.

## Elementary Covers

```python
from sage.all import GF, PolynomialRing
from ascovers import SuperellipticCurve, elementary_cover

F = GF(5)
R = PolynomialRing(F, "x")
x = R.gen()

C = SuperellipticCurve(x**3 + x + 1, 2, prec=60)
Y = elementary_cover([C.x], prec=60)

print(Y.group)
print(Y.rhs)
print(Y.z[0]**5 - Y.z[0] == Y.x)
```

Here `Y` is the cover with equation `z0^5 - z0 = x`.

## Functions and Forms

`ArtinSchreierFunction` represents rational expressions in the cover
coordinates `x, y, z_i`.  The defining equations are available through
`reduce()`:

```python
z = Y.z[0]
print((z**5).reduce())     # z0 + x
print(z.group_action(Y.group.gens[0]))
```

`ArtinSchreierForm` represents `h dx`.  Its coefficient is an
`ArtinSchreierFunction`, so local expansions and group actions are shared with
the function class:

```python
omega = z * Y.dx
print(omega.form)
print(omega.expansion_at_infty())
```

## Templates and Groups

The migrated package currently exposes templates for:

- elementary abelian covers: `elementary_template`, `elementary_cover`,
- Artin-Schreier-Witt covers: `witt_template`, `witt_cover`,
- Heisenberg covers: `heisenberg_template`, `heisenberg_cover`,
- quaternion covers in characteristic two: `quaternion_template`,
  `quaternion_cover`,
- dihedral/D8-style templates: `dihedral_template`, `dihedral_cover`,
  `d8_template`, `d8_cover`.

Group constructors include `cyclic_group`, `elementary_abelian_group`,
`heisenberg_group`, `quaternion_group`, and `dihedral_group`, with legacy names
such as `cyclic_gp`, `elementary_gp`, and `dihedreal_gp`.

## Cech-De Rham Cocycles

`ArtinSchreierDeRhamCocycle` represents a cocycle `(omega0, f, omega8)` with
`omega8 = omega0 - df`.

```python
cocycle = Y.cech(z.diffn(), z)
print(cocycle.omega8 == 0 * Y.dx)
```

## Polydifferentials and Group Actions

The package also includes the migrated helpers for holomorphic
polydifferentials, group-action matrices, canonical-ideal components, and
Magma command generation:

```python
holo_basis = Y.holomorphic_differentials_basis()
matrices = Y.group_action_matrices_holo(basis=holo_basis)

poly_basis = Y.holo_polydifferentials_basis(2)
poly_matrices = Y.group_action_matrices_poly(2, basis=poly_basis)
```

The Magma helpers can return command text without requiring Magma locally:

```python
from ascovers import magma_module_decomposition

command = magma_module_decomposition(matrices[0], matrices[0], text=True)
```
