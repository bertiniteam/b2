# The junk-success bug: the Cauchy endgame is blind to poles

*Found 2026-07-03 during the structured-output-directory arc, via the records' own
`function_residual` field.  Present in Bertini 1 (verified, reproduction below) and
inherited faithfully by bertini2.  Companion diagram:
`junk-success-bug.puml` / `.png` beside this file.  Bertini 1 reproduction files:
`junk-success-b1-repro-input` / `-start`.*

## Symptom

On an **unpatched affine user homotopy** toward a **deficient target** — a branch with
no finite root to reach — the Cauchy endgame returns `SuccessCode::Success` at a point
that is **not a root**: function residual ~ 1.0, condition number ~ 1e11.  The
zero-dim/HomotopySolver surface then reports it as a solution (its finite-only filter
never looks at the residual).  The powerseries endgame on the identical path truncates
correctly (`SecurityMaxNormReached`).

Plain (homogenized + patched) zero-dim solves are **not** affected: with a patch, a
deficient branch becomes a *convergent projective path* to a point at infinity, the
dehomogenized norm blows up, and the security machinery sees it.  The pathology needs
the raw affine formulation, where "nowhere to go" manifests as a Laurent pole.

## The mechanism, exactly

Take H = (1−t)·B + γt·A with

    A = {x² − 1,  x·y − 1}          (two finite roots)
    B = {x² − 1,  y·(x+1) − 1}      (ONE finite root: (1, 1/2))

On the branch x = −1 the second equation gives, in closed form,

    y(t) = −1 − (1−t)/(γt)          — a simple pole, y ~ −1/(γt) as t → 0.

The Cauchy endgame estimates the path's limit as the **mean of samples around a
circle** |t| = r.  For uniformly spaced samples t_j = r·e^(2πij/N):

    mean_j [ a/t_j + b ]  =  (a/r)·(Σ_j e^(−2πij/N))/N + b  =  b,   EXACTLY,

because the roots-of-unity sum is identically zero.  So the approximation of this
*diverging* path is the **stationary finite constant** b = −1 + 1/γ, at every radius,
in exact arithmetic.  Consequences:

1. **The acceptance criterion fires.**  Convergence = two consecutive approximations
   agree to FinalTolerance.  They are the same constant; they agree to refinement
   noise (~1e−12 observed).  The criterion itself is sound — the *estimator* violates
   the axiom behind it ("a diverging path can never produce agreeing estimates"),
   because its kernel contains exactly the divergent signal.
2. **The security check could never fire**, because it watched the norm of the
   **extrapolated approximation** (the finite mean b), not the samples.  (Fixed:
   it now watches the minimum dehomogenized norm over the loop samples.)
3. Sample-norm watching alone is still insufficient *for this path*: acceptance
   happens at the second loop (|t| ≈ 0.05, sample norms ~12–23), while max_norm (1e4)
   is first exceeded at |t| ≈ 1.2e−5 — thirteen radius-halvings after the endgame
   has returned.

Powerseries is immune because its Hermite extrapolation *follows the samples* — a pole
sends the extrapolation large, agreement fails, and security truncates.

## Reproduction, bertini2 (Python)

```python
import bertini as pb
from bertini import Variable, VariableGroup, System
from bertini.nag_algorithm import blend_homotopy

x, y = Variable('x'), Variable('y')
def sys_of(*fns):
    s = System(); s.add_variable_group(VariableGroup([x, y]))
    for f in fns: s.add_function(f)
    return s

A = sys_of(x**2 - 1, x*y - 1)
B = sys_of(x**2 - 1, y*(x + 1) - 1)      # deficient: one finite root

r1 = pb.solve(A, seed=42, directory='repro')
r2 = pb.solve(B, homotopy=blend_homotopy(B, A), start=r1, seed=42, directory='repro')
for m in r2.solver.solution_metadata():
    print(int(m.path_index), float(m.function_residual))   # bug: one path ~1.0, Success
```

Pinned tests:
- C++ (the correctness gate): `deficient_affine_user_homotopy_never_junk_success` in
  `core/test/endgames/generic_cauchy_test.hpp` — runs in all 8 endgame variants.
- Python symptom pin: `test_chained_deficient_target_paths_never_junk_success` in
  `python/test/records/records_test.py`.

## Reproduction, Bertini 1 (verified on v1.7)

`input`:

```
CONFIG
UserHomotopy: 1;
EndgameNum: 2;
SecurityLevel: 0;
TrackType: 0;
END;
INPUT
variable x, y;
function f1, f2;
pathvariable t;
parameter s;
s = t;
gamma = 4/5 + 3/10*I;
f1 = (1-s)*(x^2-1) + gamma*s*(x^2-1);
f2 = (1-s)*(y*(x+1)-1) + gamma*s*(x*y-1);
END;
```

`start`:

```
2

1 0;
1 0;

-1 0;
-1 0;
```

Run `bertini input start`.  Observed (main_data), path 1:

    endpoint          (−1,  0.0958904109... − 0.4109589041...i)
    function residual 1.0
    condition number  2e28
    accuracy estimate 1.09e−12  (consecutive endpoint estimates agreed)
    counted           "there appear to be 2 solutions", multiplicity 1

The endpoint is **exactly** −1 + 1/γ = the Laurent constant the theory predicts — B1's
Cauchy mean annihilated the pole identically.  B1 also tracked to T ≈ 1e−4, where the
true path norm exceeds its SecurityMaxNorm, without truncating: its security check
also watches the extrapolation.  With `EndgameNum: 1` (powerseries) B1 prints
"Truncated infinite paths: 1" — correct, matching bertini2's powerseries.

## Why the affine case is pernicious (and projective is fine)

With SecurityLevel 1 on a homogenized+patched system, the points at infinity are
DELIBERATELY computed: in projective space, infinity is just a finite point on the
patch, and the endgame converges to it like any other.  Truncation there is an
optimization.  In a raw AFFINE homotopy, infinity is really infinity -- the path and
the endgame can never converge, so truncation is the only honest outcome, and an
endgame that manufactures a finite "limit" (the pole's circle-mean) is manufacturing
an answer where none exists.

## A second manifestation: singular targets, misdetected cycle

Chaining into a target with a multiplicity-2 root (e.g. {(x-1)^2+(y-1)^2, x-y} from
the same family) shows a sibling failure: both paths genuinely converge to (1,1)
(verified by independent continuation; the approach is ~sqrt(t), i.e. cycle 2), but
the endgame detects CYCLE 1, its one-circuit loops are not closed on the two-sheeted
branch, and the resulting garbage means stabilize at non-roots (residuals ~3.4/1.2)
-- junk Success again.  Note for the fix: the twisted mean c_-1 is nonzero for ANY
improperly-closed loop content (not only poles), so the operating-zone check below
should catch this manifestation as well.

## Fix status

**The fix lives in its OWN PR** (user decision, 2026-07-03): discovered during the
records arc (PR #68), important enough to land independently.  PR #68 carries only
this documentation, the reproductions, and the xfailed Python pin.

- **For the fix PR -- security check watches the SAMPLES — minimum dehomogenized norm
  over the current loop, in both `RunImpl` and the AMP driver
  (`MinLoopSampleNorm` / `MinLoopSampleNormAMP`; patch staged).  Necessary
  (the old check was structurally unable to fire on a pole) but not sufficient alone:
  acceptance happens at radii where the pole's sample norms are still tiny.
- **Rejected**: a residual gate on acceptance (`ApproximationIsVerifiedRoot`,
  briefly on the branch, reverted).  Function values are scaling-sensitive — a system
  scaled by 1e6 would make true solutions fail the gate.  Endgame convergence is
  deliberately not conditional on function values.
- **Proposed, pending decision**: the pole-component operating-zone check.  The same
  loop samples give the t^(−1) Laurent coefficient by one twisted mean:

      c₋₁ = mean_j [ x_j · e^(+iθ_j) ]

  (every analytic term — constant, integer and fractional positive powers over the
  closed c-circuit loop — contributes zero, by the same orthogonality that hides the
  pole from the plain mean).  A significant ‖c₋₁‖ relative to
  FinalTolerance·(1 + ‖mean‖) means the loop contains a pole ⇒ **not in the endgame
  operating zone** ⇒ do not accept, keep advancing ⇒ the (now sample-watching)
  security check truncates at max_norm.  Coordinates only; acceptance criterion
  untouched; a converging cycle-c path has c₋₁ → 0 like r^(1/c).

## Worth reporting upstream

The B1 reproduction is the two files above; the maintainers may want it even if only
as a documentation caveat for `UserHomotopy: 1` + `EndgameNum: 2`.
