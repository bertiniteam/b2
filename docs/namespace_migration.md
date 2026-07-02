# Migrating to the flattened `bertini` API (3.0)

Ahead of the public release we flattened the Python API so the everyday things live at the top
level (`bertini.something`) and are tab-completable, while the rare/internal machinery stays nested.
This is a **breaking** change to spellings — the objects and behavior are the same, only the import
paths and a few names changed.

> Status: DRAFT skeleton. Spellings below reflect the landed renames; a final pass will confirm each
> and add examples before this goes to collaborators.

## The short version

- Common classes are now at the top level: `bertini.ZeroDimSolver`, `bertini.AMPTracker`,
  `bertini.Slice`, `bertini.SuccessCode`, ...
- The symbolic system is `bertini.symbolics` (was `bertini.function_tree`), and it is **flat** —
  `bertini.symbolics.Variable`, not `bertini.symbolics.symbol.Variable`.
- Multiprecision numbers are lowercase `_mp` (like numpy dtypes): `bertini.complex_mp`,
  `bertini.real_mp`, `bertini.int_mp`, `bertini.rational_mp`.
- `bertini.linalg` is gone: its system-building calls are now **methods on `System`**, and its
  coefficient helpers are top-level `bertini.coefficient` / `bertini.coefficients`.

## Renames (old → new)

### Symbolic system: `function_tree` → `symbolics` (now flat)
| old | new |
|---|---|
| `bertini.function_tree` | `bertini.symbolics` |
| `bertini.function_tree.symbol.Variable` | `bertini.symbolics.Variable` (or `bertini.Variable`) |
| `bertini.function_tree.symbol.Complex` / `Integer` / `Rational` | `bertini.symbolics.Complex` / `Integer` / `Rational` |
| `bertini.function_tree.symbol.E` / `Pi` | `bertini.symbolics.E` / `Pi` (or constants `bertini.E` / `bertini.Pi` / `bertini.I`) |
| `bertini.function_tree.root.NamedExpression` | `bertini.symbolics.NamedExpression` (or `bertini.Named`) |
| `bertini.function_tree.operator.Sum` / `Mult` / ... | `bertini.symbolics.Sum` / `Mult` / ... |

### Multiprecision numbers → lowercase `_mp` (top level + in `multiprec`)
| old | new |
|---|---|
| `bertini.multiprec.Complex` | `bertini.complex_mp` (also `bertini.multiprec.complex_mp`) |
| `bertini.multiprec.Float` | `bertini.real_mp` |
| `bertini.multiprec.Int` | `bertini.int_mp` |
| `bertini.multiprec.Rational` | `bertini.rational_mp` |

(`Float` was a real number — `real_mp` says so. The symbolic `Complex`/`Integer`/`Rational` keep
CapWords in `symbolics`; lowercase = a concrete number, CapWords = a symbol.)

### Hoisted to the top level (still in their submodules too)
| old | new |
|---|---|
| `bertini.nag_algorithm.ZeroDimSolver` | `bertini.ZeroDimSolver` |
| `bertini.nag_algorithm.HomotopySolver` | `bertini.HomotopySolver` |
| `bertini.nag_algorithm.SolutionPathCollector` | `bertini.SolutionPathCollector` |
| `bertini.nag_algorithm.Slice` | `bertini.Slice` |
| `bertini.tracking.AMPTracker` / `DoublePrecisionTracker` / `MultiplePrecisionTracker` | `bertini.AMPTracker` / ... |
| `bertini.tracking.SuccessCode` / `Predictor` | `bertini.SuccessCode` / `bertini.Predictor` |
| `bertini.nag_algorithm.StartSystemType` | `bertini.StartSystemType` |
| (enums generally) | at the root: `SuccessCode`, `Predictor`, `MonomialOrder`, `StartSystemType` |

### `linalg` dissolved → `System` methods + top-level helpers
| old | new |
|---|---|
| `linalg.add_functions(sys, exprs)` | `sys.add_functions(exprs)` |
| `linalg.add_linear(sys, A, x, b)` | `sys.add_linear(A, x, b)` |
| `linalg.add_linear_forms(sys, C)` | `sys.add_linear_forms(C)` |
| `linalg.add_products_of_linears(sys, F)` | `sys.add_products_of_linears(F)` |
| `linalg.add_slices_as_products(sys, slices)` | `sys.add_slices_as_products(slices)` |
| `linalg.randomize(sys, R)` | `sys.randomize(R)` |
| `linalg.coefficient(v)` | `bertini.coefficient(v)` |
| `linalg.as_coefficients(A)` | `bertini.coefficients(A)` (dropped the `as_`) |
| `linalg.slice_from_coefficients(coeffs, vars)` | `bertini.Slice.from_coefficients(coeffs, vars)` |
| `linalg.variable_vector('v', n)` | `bertini.variables('v', n)` (list) or `bertini.VariableGroup('v', n)` |

### Endgame class names spelled out (`EG` → `Endgame`, `PS` → `PowerSeries`)
| old | new |
|---|---|
| `endgame.AMPCauchyEG` | `endgame.AMPCauchyEndgame` |
| `endgame.AMPPSEG` | `endgame.AMPPowerSeriesEndgame` |
| `endgame.FixedDoubleCauchyEG` / `FixedDoublePSEG` | `endgame.FixedDoubleCauchyEndgame` / `FixedDoublePowerSeriesEndgame` |
| `endgame.FixedMultipleCauchyEG` / `FixedMultiplePSEG` | `endgame.FixedMultipleCauchyEndgame` / `FixedMultiplePowerSeriesEndgame` |

### Settings
- Set individual settings by name in one call: `owner.set(**kwargs)` (an alias of the existing
  `owner.update(**kwargs)`), e.g. `solver.set(final_tolerance="1e-11")`. Config classes stay in
  their submodules (`tracking.SteppingConfig`, `endgame.SecurityConfig`, ...).

## Removed — use instead
| removed | use instead |
|---|---|
| `bertini.container` (module) | `bertini.VariableGroup` (top level) |
| `bertini.multiprec.Vector` | a plain numpy array |
| `linalg.variable_matrix` | build a numpy object array of `Variable`s yourself |
| `function_tree.symbol.make_e` / `make_i` / `make_pi` | the constants `bertini.E` / `bertini.I` / `bertini.Pi` |
| `symbolics.Differential` (hidden) | differentiate with an explicit variable: `f.differentiate(x)` / `bertini.jacobian(...)` |
| NID classes + `WitnessSetMultiplePrecision` (hidden) | not yet public — NID's `Solve()` is not implemented |

## Notes
- During the transition several old spellings still resolve as **deprecated aliases**
  (`bertini.function_tree`, `bertini.multiprec.Complex`, `endgame.AMPCauchyEG`, ...); they are
  unadvertised and will be removed. Migrate to the new spellings.
- The auto-generated API reference re-walks the package, so it already reflects the flat surface.
