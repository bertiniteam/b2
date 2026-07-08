# ADR-0050: Publish the versioned docs by serving the `docs-store` branch, not `actions/deploy-pages`

**Status:** Accepted
**Date:** 2026-07-08

## Context

bertini2.org keeps historical docs: each real release lives forever under `/vX.Y.Z/`, a
moving `/stable/` points at the newest release, and the root is a landing page listing
versions (see the docs-versioning scheme). `build_docs.yml` builds a single-version `site/`
(Doxygen → `site/cpp`, Sphinx **furo** → `site/python`, custom landing → `site/index.html`)
and `tools/assemble_versioned_docs.py` folds it into a persistent **`docs-store` branch**:
the durable byte store, one `/vX.Y.Z/` per release, `versions.json` + root `index.html`
derived from whatever version dirs are present (truth-vs-derived, mirroring the records
doctrine).

The original plan (and the first implementation) then published that store to GitHub Pages
with `actions/upload-pages-artifact` + `actions/deploy-pages` — Pages `build_type: workflow`
— deliberately to avoid a Pages *settings* change. **This did not work for the multi-version
store**, and the failure was expensive to diagnose:

- The first bug was a gate: `deploy_docs` keyed on a cross-job output
  (`needs.build_docs.outputs.is_release`) that **did not propagate** into the downstream job
  `if` in the reusable-workflow context. The job *skipped* even though `build_docs` had set
  the output and uploaded the artifact.
- Fixing the gate exposed the real wall: `actions/deploy-pages` then failed **structurally**
  with an Azure **`BlobNotFound`** ("the specified blob does not exist"), *reproducibly*,
  across re-runs. Yet:
  - the `github-pages` artifact was valid and **downloaded fine** via the normal artifacts
    API (a well-formed ~140 MB / 9.6k-file site: `.nojekyll`, root `index.html`, no
    symlinks);
  - permissions were not the cause (a flat single-version deploy in **May 2026 succeeded
    with fewer permissions**);
  - it was not transient (identical instant failure on every retry).

The one thing that differed from the working May deploy was the **payload** — a flat
single-version site (worked) vs. the ~140 MB versioned store (failed). Whatever the precise
trigger inside GitHub's Actions→Pages *artifact-exchange*, it is opaque to us, only
reproducible through ~20-minute CI cycles, and not something we control or can pin.

## Decision

**Stop deploying Pages from an Actions artifact. Serve GitHub Pages directly from the
`docs-store` branch** (Pages `build_type: legacy`, source = `docs-store` `/`).

Because `build_docs` already **pushes** the assembled store to `docs-store`, that push *is*
the deploy — GitHub rebuilds the branch-source Pages site automatically. Consequences:

- **`actions/deploy-pages` and `actions/upload-pages-artifact` are removed**, and the entire
  `deploy_docs` job is deleted. There is no artifact exchange, so `BlobNotFound` cannot
  occur; there is no cross-job output feeding a downstream `if`, so the propagation bug is
  moot. `publish.yml`'s docs caller drops the now-unused `pages:` / `id-token:` permissions
  and keeps only `contents: write` (to push the branch).
- **The custom domain moves into the branch.** Switching the Pages source *cleared* the
  `bertini2.org` custom domain, because branch-source Pages reads the domain from a `CNAME`
  file. We add `/CNAME` to `docs-store` and re-set the Pages `cname`; the assembler writes
  `/CNAME` via `--cname bertini2.org` so every rebuild preserves it. **This ADR knowingly
  overrides the earlier "no Pages settings change" constraint** — that constraint assumed
  the Actions deploy worked, and it did not.
- **The published store is slimmed** (independently useful, and it keeps the branch small
  since it grows per release): `/stable/` becomes a small **redirect stub** to the newest
  `/vX.Y.Z/` instead of a 69 MB byte-for-byte copy, and Sphinx `.doctrees` caches are
  dropped (`sphinx-build -d` puts them outside `site/`, and the assembler ignores them on
  copy as defense-in-depth). Store: 173 MB → 86 MB, 9589 → 4597 files.

## Consequences

- **Docs updates are decoupled from releases.** A prose/source fix = dispatch `build_docs`
  (it assembles into `docs-store` and pushes → Pages republishes); a built-HTML typo or a
  rollback = push `docs-store` directly. No PyPI, no wheel matrix, no `deploy-pages`.
- **Do not switch docs publishing back to `actions/deploy-pages`.** It failed structurally
  with `BlobNotFound` on the versioned store and is a GitHub-side black box; branch-source
  is the working path. If you must revisit it, keep branch-source live until a full
  multi-version deploy is *proven* green.
- **`/stable/` deep links bounce.** A link to `/stable/some/page` redirects to the version
  root, not the exact page — the right trade for a "current docs" entry point; the landing
  card links straight to `/vX.Y.Z/`.
- **Tooling/source coupling caveat.** `build_docs` checks out `inputs.ref` for *both* the
  doc source and the tooling, so re-publishing an *already-released* version's store with
  *newer* tooling can't go through the workflow (the old tag carries the old assembler —
  e.g. no `--cname`). Re-slim such a store by hand for one-offs; the pipeline runs clean for
  future releases, whose tags carry matching tooling.
- **The custom domain now depends on the `/CNAME` file** in `docs-store`. Deleting it (or
  regenerating the branch without `--cname`) drops `bertini2.org` on the next Pages build.
