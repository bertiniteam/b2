# Migrating to the flattened `bertini` API (3.0)

Ahead of the public release we flattened the Python API so the everyday things live at the top
level (`bertini.something`) and are tab-completable, while the rare/internal machinery stays nested.
This is a **breaking** change to spellings — the objects and behavior are the same, only the import
paths and a few names changed. The old spellings have been **removed** (no deprecation window), so
update call sites in one pass; the table below is the map.

## The short version

- Common classes are now at the top level: `bertini.ZeroDimSolver`, `bertini.AMPTracker`,
  `bertini.Slice`, `bertini.SuccessCode`, ...
- The symbolic system is `bertini.symbolics` (was `bertini.function_tree`), and it is **flat** —
  `bertini.symbolics.Variable`, not `bertini.function_tree.symbol.Variable`.
- Multiprecision numbers are lowercase `_mp` (like numpy dtypes): `bertini.complex_mp`,
  `bertini.real_mp`, `bertini.int_mp`, `bertini.rational_mp`.
- `bertini.linalg` is gone: its system-building calls are now **methods on `System`**, and its
  coefficient helpers are top-level `bertini.coefficient` / `bertini.coefficients`.

## Before / after

```python
# --- BEFORE ---
import bertini as pb
from bertini import linalg
from bertini.nag_algorithm import ZeroDimSolver
from bertini.multiprec import Complex

x = linalg.variable_vector('x', 3)              # numpy array of variables
lam = pb.function_tree.symbol.Variable('lam')
A = np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]])

sys = pb.System()
sys.add_variable_group(pb.container.VariableGroup([*x, lam]))
linalg.add_functions(sys, A @ x - lam * x)      # 3 eigen-equations
sys.add_function(x[0] - 1)                       # normalize (square: 4 eqns, 4 vars)
gamma = Complex('0.6', '0.8')

solver = ZeroDimSolver(sys)
solver.solve()
```

```python
# --- AFTER ---
import bertini as pb

x = np.array(pb.variables('x', 3), dtype=object)   # variables() returns a list; wrap for `@`
lam = pb.Variable('lam')
A = np.array([[2, 1, 0], [1, 3, 1], [0, 1, 4]])

sys = pb.System()
sys.add_variable_group(pb.VariableGroup([*x, lam]))
sys.add_functions(A @ x - lam * x)                 # a System method now
gamma = pb.complex_mp('0.6', '0.8')

solver = pb.ZeroDimSolver(sys)                     # hoisted to the top level
solver.solve()
```

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

The `.symbol` / `.root` / `.operator` sub-levels are gone — everything is directly in `symbolics`.

### Multiprecision numbers → lowercase `_mp` (top level + in `multiprec`)
| old | new |
|---|---|
| `bertini.multiprec.Complex` | `bertini.complex_mp` (also `bertini.multiprec.complex_mp`) |
| `bertini.multiprec.Float` | `bertini.real_mp` |
| `bertini.multiprec.Int` | `bertini.int_mp` |
| `bertini.multiprec.Rational` | `bertini.rational_mp` |

`Float` was a real number — `real_mp` says so. Casing is a signal: **lowercase** = a concrete number,
**CapWords** in `symbolics` = a symbol (the symbolic `Complex`/`Integer`/`Rational` keep CapWords).

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

The homotopy helpers (`parameter_sweep`, `moving_homotopy`, `coefficient_parameter_homotopy`,
`blend_homotopy`) and the `*Config` structs stay in `nag_algorithm` / `tracking` / `endgame`.

### `linalg` dissolved → `System` methods + top-level helpers
| old | new |
|---|---|
| `linalg.add_functions(sys, exprs)` | `sys.add_functions(exprs)` |
| `linalg.add_linear(sys, A, x, b)` | `sys.add_linear(A, x, b)` |
| `linalg.add_linear_forms(sys, C)` | `sys.add_linear_forms(C)` |
| `linalg.add_products_of_linears(sys, F)` | `sys.add_products_of_linears(F)` |
| `linalg.add_slices_as_products(sys, slices)` | `sys.add_slices_as_products(slices)` |
| `linalg.randomize(sys, R)` | `sys.randomize(R)` (or `sys.randomize()` for a generic R) |
| `linalg.coefficient(v)` | `bertini.coefficient(v)` |
| `linalg.as_coefficients(A)` | `bertini.coefficients(A)` (dropped the `as_`) |
| `linalg.slice_from_coefficients(coeffs, vars)` | `bertini.Slice.from_coefficients(coeffs, vars)` |
| `linalg.variable_vector('v', n)` | `bertini.variables('v', n)` — see note |

**`variable_vector` note:** the old helper returned a **numpy object array**; `bertini.variables`
returns a **list**. If you feed it to numpy (`A @ x`, `lam * x`), wrap it:
`np.array(bertini.variables('v', n), dtype=object)`. For a variable *group*,
`bertini.VariableGroup('v', n)` builds `v0..v{n-1}` directly.

`numpy.linalg` is unaffected — `np.linalg.norm` / `np.linalg.eigvals` were never `bertini.linalg`.

### Endgame class names spelled out (`EG` → `Endgame`, `PS` → `PowerSeries`)
| old | new |
|---|---|
| `endgame.AMPCauchyEG` | `endgame.AMPCauchyEndgame` |
| `endgame.AMPPSEG` | `endgame.AMPPowerSeriesEndgame` |
| `endgame.FixedDoubleCauchyEG` / `FixedDoublePSEG` | `endgame.FixedDoubleCauchyEndgame` / `FixedDoublePowerSeriesEndgame` |
| `endgame.FixedMultipleCauchyEG` / `FixedMultiplePSEG` | `endgame.FixedMultipleCauchyEndgame` / `FixedMultiplePowerSeriesEndgame` |

### Settings
Set individual settings by name in one call: `owner.set(**kwargs)` — an alias of the existing
`owner.update(**kwargs)`, e.g. `solver.set(final_tolerance="1e-11")` or
`tracker.get_stepping().set(max_step_size="0.05")`. The `*Config` classes stay in their submodules
(`tracking.SteppingConfig`, `endgame.SecurityConfig`, `nag_algorithm.TolerancesConfig`, ...).

## Removed — use instead
| removed | use instead |
|---|---|
| `bertini.container` (module) | `bertini.VariableGroup` (top level) |
| `linalg.variable_matrix` | a numpy object array you build (e.g. `np.array([[Variable(f'm_{i}_{j}') ...]], dtype=object)`) |
| `function_tree.symbol.make_e` / `make_i` / `make_pi` | the constants `bertini.E` / `bertini.I` / `bertini.Pi` |
| `symbolics.Differential` (hidden) | differentiate with an explicit variable: `f.differentiate(x)` / `bertini.jacobian(...)` |
| NID classes + `WitnessSetMultiplePrecision` (hidden) | not yet public — NID's `Solve()` is not implemented |

`bertini.multiprec.Vector(n)` is **kept** (a pre-sized numpy array of the multiprecision-complex
dtype, handy for result buffers); `np.zeros(n, dtype=bertini.complex_mp)` is the explicit equivalent.

## Notes
- The auto-generated API reference re-walks the package, so it already reflects the flat surface; the
  new hand-written **Everyday API** page groups the common names by task.
- If an old spelling still appears somewhere, it will now raise `AttributeError` / `ImportError` —
  the map above has its replacement.
