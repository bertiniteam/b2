# Ledger v0 — pure-Python prototype

A walk-away-able pilot of the solution provenance ledger
(`arcs/solution-provenance-ledger.md`): plain-text records, memoized solving, and
invisible resume — built entirely on the public Python API (no C++ changes beyond what
PR #64 already shipped, whose `content_digest()` this leans on).

**This directory is deliberately outside `python/bertini/`** so nothing here ships in the
wheel.  Abandoning the experiment is `rm -rf prototypes/ledger_v0` + deleting the branch.

## What it demonstrates

- **Object store**: content-addressed files under `objects/<2hex>/<digest>` — definitions
  (target system classic input, pickled homotopy instances) stored once, referenced by id.
- **Journals**: per-run append-only JSONL under `journals/` — one record per line,
  readable with grep/jq/pandas, torn-last-line tolerant.
- **`ensure_solved(target, ledger)`** — solve() as "ensure answered":
  - a completed run for this ask → returns recorded solutions (no tracking);
  - a partial run (crashed) → adopts the recorded homotopy instance and finishes only
    the missing paths;
  - nothing → creates the run and computes everything.
- **The kill-and-rerun pilot**: `demo_resume.py` solves with a simulated walltime kill
  partway, then reruns and finishes from where it died.

## v0 stopgaps (known, deliberate)

- **Instance persistence by pickle, not by seed.**  Seed-rooted randomness derivation
  (arc rung 2) does not exist yet, so the homotopy instance (gamma, start-system
  coefficients) is persisted as a pickled blob in the object store and re-adopted on
  resume — the "run manifest" mechanism.  Once derivation lands, the blob becomes
  redundant and the ask key gains `seed`.
- Pickle blobs are binary — the one non-text artifact.  They vanish with seed derivation.
- Object files are named by `System.content_digest()` but contain the *classic input*
  rendering (the canonical encoding text is not yet exposed to Python), so the
  self-verifying property (hash the file, get its name) does not hold in v0.
- Single-process only: no manager/worker funnel, no fsync batching.  The journal format
  and one-writer-per-journal rule are the HPC-relevant parts being piloted.

## Record schema (`ledgerrec/0`)

Every line is a JSON object with a `kind`:

- `{"kind": "run", "schema": "ledgerrec/0", "run": <run id>, "ask": {"target": <digest>,
   "config": {...}}, "homotopy_blob": <object id>, "target_object": <object id>,
   "num_paths": N, "start_points": [[[re,im],...], ...]}` — the start points are
   recorded as text in the header (they ARE provenance: the canonical start labels'
   values), so resume needs no start-system object.
- `{"kind": "track", "run": <run id>, "index": i, "status": "success"|"failed",
   "endpoint": [[re_str, im_str], ...] | null}`

## Run it

```bash
PYTHONPATH=python pytest prototypes/ledger_v0/test_ledger.py
PYTHONPATH=python python prototypes/ledger_v0/demo_resume.py
```
