# ADR-0061: at most one encoding version per released library version

**Status:** Accepted
**Date:** 2026-09-19

## Context

Systems and configs have stable cross-session identities built on versioned canonical
encodings, `b2sysenc/<n>` (ADR-0042) and `b2cfgenc/<n>` (ADR-0043).  The version token is
part of the digest preimage, so moving it makes every previously written record a different
ask.  A registry file records every version there has ever been, and a test refuses a
regenerated fixture that skipped the bump.

The rule for the registry was "append-only, never edit a line".  That is the right instinct
about *shipped* history and the wrong rule for development, because it makes the version
climb with commits rather than with releases.  The evidence is in the file: `b2cfgenc/1` and
`/2` were both minted while 3.0.0 was being developed and **neither ever shipped** — v3.0.0,
the first release, already wrote `/3`.  So the number said three encodings existed when one
had ever been used, and nothing in the file said which.

Two encoding changes are queued for 3.5 (#402 and #403).  Under append-only, an unreleased
3.5 would go `/3` → `/6`, and only `/6` would ever exist in the world.

A second, separate defect: the registry hashes the *keyspace* — what each config encoder
turns one recipe into.  It cannot see the **composition**, which configs a solve folds into
its settings text.  `b2cfgenc/3` and `/4` carry byte-identical keyspace hashes, because `/4`
*was* MidPathConfig joining the ask and no encoder changed.  Those two lines record "nothing
changed" about a bump that moved every settings digest in existence.

## Decision

**At most one encoding version per released library version.**

- Every registry line records the release it first shipped in, or `never-shipped` for one
  superseded before any release used it.  The two historical never-shipped `b2cfgenc` lines
  are labelled as such rather than quietly implying they were real.
- While a release is **unreleased**, an identity-affecting change updates the top registry
  line **in place** — new hash, same version token — and regenerates the golden fixture.
  Once a version has shipped, its line is history and is never edited: the bump-and-append
  path applies from then on.  The failing test works out which case applies and says so.
- No two lines may claim the same release, and releases increase down the file.  A pre-release
  suffix is not a release of its own: `3.5.0.dev1` and `3.5.0` are the same release series, so
  a dev wheel cannot buy another version token.
- The `b2cfgenc` registry gains a **composition** column: a hash of the ordered list of configs
  a zero-dim solve folds into its settings text.  It is owned by the test that can build a
  solver (`zero_dim_records`), while the keyspace column stays with the test that cannot
  (`config_digest_test`) — one file, two columns, two owners.

The rule and its messages live in `core/test/utility/encoding_registry.hpp`, as a pure
function over registry lines so the rule itself is testable (`encoding_registry_test.cpp`)
rather than only exercised against files that are supposed to be correct.

## Consequences

- The version number means "the *n*th encoding b2 has ever written", and a reader holding a
  record can name the release that wrote it.  Do **not** restore plain append-only: it makes
  the number count development events, which is what produced two phantom versions already.
- Records written by an earlier build of the *same unreleased* release line stop being
  recalled when its line is edited.  This is deliberate and it is cheap: the digest changes,
  so the run is a new ask and is recomputed, and nothing is corrupted or silently reused.  A
  development line behaves like a development line.
- Editing a line is only legitimate while it is unreleased.  The test cannot see the previous
  contents of the file, so this last part is discipline backed by review, not by machinery —
  which is why the release column exists to make the claim visible in the diff.
- The composition column means a composition-only bump now leaves evidence.  Do not remove it
  on the grounds that the keyspace hash already covers the encoders; it demonstrably does not
  cover this, which is how MidPathConfig stayed outside the ask.
- **Renumbering the release under development means editing the release column of every line
  that has not shipped, and nothing enforces it.**  A last line naming an earlier release is a
  perfectly ordinary state — it is what the file looks like whenever a release has shipped and
  nothing has changed the encoding since — so the machinery cannot tell that case apart from a
  stale column left behind by a renumber.  Measured when 3.5 became 4.0: with the column still
  reading `3.5.0`, every test passes, and the *next* encoder change is then told to append a new
  version token instead of editing the line in place, minting a version for an encoding that
  never shipped.  That is the exact mistake this ADR exists to prevent, arriving by a different
  door.  Renumber the columns in the same commit as the `VERSION` file.
