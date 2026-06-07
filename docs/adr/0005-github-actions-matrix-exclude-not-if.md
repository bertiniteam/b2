# ADR-0005: Use matrix exclude to filter CI jobs, not matrix.os in job-level if

**Status:** Accepted  
**Date:** 2026-06-07

## Context

The `test_wheels_linux_macos` job runs the full Python test suite on host runners
after the wheel builds. Linux wheels are tested inside the manylinux container
(see ADR-0003), so this job should only run for macOS — ubuntu host runner tests
would crash on the manylinux wheel.

The first implementation used a job-level `if:` expression:

```yaml
test_wheels_linux_macos:
  if: ${{ matrix.os != 'ubuntu-latest' && inputs.minimal != true }}
```

This caused an immediate 0-second "workflow file issue" failure. All job names in
the run appeared as unexpanded expressions (`${{ matrix.os }} Python-...` instead
of the actual OS name), which is GitHub Actions' symptom for a workflow validation
error before execution begins.

### Why the job-level if fails

GitHub Actions validates job-level `if:` expressions **before** the matrix is
expanded. When `matrix.os` is dynamically generated via:

```yaml
matrix:
  os: ${{ fromJson(needs.set_matrix.outputs.os) }}
```

GitHub has no information about what values `matrix.os` can take at validation time.
The expression `matrix.os != 'ubuntu-latest'` therefore fails validation — GitHub
cannot determine whether it is syntactically valid to reference `matrix.os` in a
job-level condition before the matrix is known.

This is a documented GitHub Actions limitation: **job-level `if:` cannot reference
`matrix.*` when the matrix is dynamically generated from a prior job's output**.
Step-level `if:` inside a job works fine (the matrix is materialized by then), but
job-level `if:` is evaluated at workflow parse time.

## Decision

Use matrix `exclude:` instead of a job-level `if:`:

```yaml
test_wheels_linux_macos:
  if: ${{ inputs.minimal != true }}    # only non-matrix conditions here
  strategy:
    matrix:
      os: ${{ fromJson(needs.set_matrix.outputs.os) }}
      python-version: ${{ fromJson(needs.set_matrix.outputs.python_versions) }}
      exclude:
        - os: ubuntu-latest
```

`exclude:` is evaluated after the matrix is materialized — GitHub simply removes
the matching combinations from the generated job list. It works correctly with
`fromJson()` matrices.

## Consequences

- **Validation passes.** GitHub can parse and validate the workflow before any runner
  starts.
- **ubuntu-latest excluded correctly.** No Python test jobs run on ubuntu host runners
  for this job. Linux wheels are verified only by the smoke test inside manylinux
  (ADR-0003).
- **Rule for future workflow authors:** Never reference `matrix.*` in a job-level
  `if:` when the matrix comes from `fromJson()` of a prior job's output. Use
  `exclude:` for OS/version filtering and reserve job-level `if:` for non-matrix
  conditions (`inputs.*`, `github.event_name`, etc.).
- **Step-level `if:` is fine.** If per-step filtering by `matrix.os` is needed, put
  the condition on the step, not the job.
