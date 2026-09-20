# ADR-0063: a setting is a config field, and derived state is rederived

**Status:** Accepted
**Date:** 2026-09-20

## Context

Every config owner -- tracker, endgame, solver -- holds its settings in `detail::Configured`
structs, reachable uniformly by `Get<T>()`/`Set<T>()`, by the flat `update(field=...)` router,
and as a bundle through `get_settings`/`set_settings`.  The settings digest that identifies a
solve is the canonical text of those structs, and the configuration reference page is generated
from them.

Four settings sat outside that.  Three were bare members of the tracker -- the predictor, the
tracking tolerance, the path truncation threshold -- each with a hand-written setter and getter.
The fourth, `final_tolerance`, was the opposite problem: it was a field of two structs at once,
the solver's `TolerancesConfig` and the endgame's `EndgameConfig`, with the solver's copy pushed
onto the endgame at solve time.

Everything built on the config surface therefore had a hole in it.  A settings bundle carried
neither the tracking tolerance nor the predictor; the settings digest did not see them, so two
solves differing only in predictor were the same ask; the generated reference could not list
them, because no amount of introspecting config structs finds a method.  `path_truncation_threshold`
had no route through a solver at all: the only way to change it was to fetch the tracker and call
a method, which is how it went unnoticed that the solver could not set it.  And a name owned by
two structs forced the flat router to pick one, making the two values agree by a push rather than
by there being one of them.

There was a reason the bespoke setters existed: each one also updated state derived from the
setting.  `SetTrackingTolerance` recomputed the digit count implied by the tolerance;
`SetPredictor` rebuilt the predictor object and re-read its order.  `Configured::Set<T>` is a
dumb store -- it overwrites the struct and tells nobody -- so a setting moved into a config and
left with that design would take effect through one route and not another.

## Decision

**A setting is a config field.**  The three tracker settings became `tracking::TrackerConfig`,
added to every tracker's `NeededConfigs`.  `final_tolerance` belongs to `EndgameConfig` alone,
the endgame being what achieves it; `TolerancesConfig` keeps the two Newton tolerances, which are
the solver's to drive.  No field name is shared by two configs.

**Derived state is rederived where it is used, not cached where it is set.**  The tracker syncs
what it derives at the top of `TrackPath` -- rebuilding the predictor object only when the choice
actually changed -- and computes the tolerance's digit count on demand.  Laziness is correctness
here: a value that arrives through `Set<TrackerConfig>`, through a restored settings bundle, or
through the flat router is as effective as one that arrives through a setter.

The named setters and getters stay, and read and write the config.  They are the terse spelling
of a common operation, not a second source of truth.

## Consequences

- Do **not** re-add a bare settings member to a config owner.  If a setting has no config, it is
  invisible to the digest, to settings bundles, and to the generated reference -- all three
  silently.
- Do **not** cache derived state in a setter as an optimization.  The config can be replaced
  wholesale by three other routes, none of which will call it.  The rebuild is guarded on the
  value having changed, so the repeated path costs a comparison.
- Do **not** give `TolerancesConfig` a `final_tolerance` back, and do not reintroduce a push from
  the solver onto the endgame.  A solver-level `final_tolerance=` still works; it routes to the
  endgame's config.
- Both moves change the settings text, so `b2cfgenc/4`'s registry line was edited in place under
  ADR-0061 -- 4.0.0 is unreleased.  The keyspace hash and the golden digest fixture moved with it.
- Four tracker settings remain method-only (`infinite_truncation`, `precision_preservation`,
  `reinitialize_initial_step_size`, and the wall-clock deadline of ADR-0060).  The first three are
  candidates for the same treatment; the deadline is deliberately out of the ask.  They are
  hand-listed in the configuration reference.
- The partition of settings across configs is still uneven -- `NewtonConfig` and `TrackerConfig`
  both describe how the tracker takes a step.  Repartitioning is deferred; this decision is about
  where a setting lives at all, not about which struct is the right home.
