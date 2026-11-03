# ADR-0074: Integration branches are first-class in CI

**Status:** Accepted
**Date:** 2026-11-03

## Context

Some lines of work should not reach `develop` until they have earned it.  The first is interval evaluation for certification (#412): it brings FLINT into the build, a dependency every source build would then need and every wheel would carry, and silviana's position is that running an experiment and developing it should not require committing the library to its dependency ("i just don't think i should have to mainline commit to a library, to run the experiment and try to develop it").

Such work needs cross-platform CI throughout, not at the end.  But CI was shaped around `develop`: the push trigger named `develop` and `main`, and `set_matrix` gave the full OS and Python matrix only to PRs based on them, one platform to everything else.  A PR into any other branch would have been verified on Ubuntu x86_64 with one Python, and merges into that branch would have run nothing at all, so the branch would never have saved the caches its PRs restore: a PR restores caches from its base branch, never from another PR, so a dependency `develop` has never built compiles cold on every PR.

The alternatives considered, and what the larger projects do: a rule naming one branch (what people do by accident, and what silviana was wary of); a per-PR label asking for the full matrix (Rust's `bors try`, Kubernetes' `/test all`), which cannot cover merges into the branch and so cannot warm caches; full CI on every PR, which only slows the small ones and still needs the push rule; and trunk-based development with the dependency behind a build option (OpenCV's `WITH_*`, PyTorch's `USE_*`), which needs the code to exist before there is anything to put behind the flag.

## Decision

A branch named `integration/<topic>` is an integration branch, and CI treats it like `develop`: the push trigger includes `integration/**`, so every merge into it runs the full matrix and saves the caches; and `set_matrix` gives a PR whose base is an integration branch the full matrix.  The rule names a kind of branch, as `release/x.y` does, not a branch.

Work lands on an integration branch through PRs from the fork, as it would on `develop`.  The branch lives on `bertiniteam/b2`, because a PR's base must be in the repository the PR is opened in.  Build dependencies the branch needs and `develop` does not (a deps image, macOS tarballs) are published under their own content-keyed tags from the branch; `develop`'s stay as they are.

Before such work reaches `develop`, a heavy new dependency becomes a build option, default off, with CI building both ways, so `develop` can carry the code without committing to the dependency; that is the trunk-based pattern, and it becomes possible once the code exists.  The first integration branch is `integration/intervals`.

## Consequences

- Do not replace `integration/**` with the name of a particular branch, and do not add a branch name to the push trigger or the matrix condition for one line of work; name the branch `integration/<topic>` instead.
- Do not take an integration branch's PRs back to the one-platform subset: the branch is held to `develop`'s bar throughout, which is the point.
- Do not "simplify" the push trigger by removing `integration/**`; without it the branch saves no caches and its every PR compiles cold.
- The Windows jobs are gated on the manual inputs, not on `set_matrix`, so they already ran on every PR; this ADR changes nothing about them.
- A `full-ci` label for one-off cross-platform runs of an ordinary branch would compose with this and is not yet added.
