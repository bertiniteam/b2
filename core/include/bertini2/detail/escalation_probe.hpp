// TEMPORARY instrumentation for the AMP precision-escalation investigation
// (branch perf/amp-block-precision-escalation).  Per-cause atomic counters so a
// solve can be attributed to WHERE precision escalations originate.  Remove (or
// gate) before folding into the PR.  See .claude/plans/typed-wondering-sutherland.md.
#pragma once

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cmath>

namespace bertini {
namespace probe {

inline bool trace_enabled()
{
	static const bool e = (std::getenv("BERTINI_TRACE_PREC") != nullptr);
	return e;
}

inline void trace_precision_change(const char* where, unsigned from, unsigned to)
{
	if (trace_enabled())
	{
		std::fprintf(stderr, "  [prec] %-18s %u -> %u\n", where, from, to);
		std::fflush(stderr); // unbuffered so a grinding process's trace survives a kill
	}
}

// thread-safe because endgames run on std::thread workers.
inline std::atomic<std::int64_t> tracker_precision_increases{0}; ///< ChangePrecision with new > current
inline std::atomic<std::int64_t> endgame_refine_escalations{0};  ///< amp_endgame RefineSampleImpl escalation branch
inline std::atomic<std::int64_t> corrector_track_hpn{0};         ///< HigherPrecisionNecessary from the tracking corrector
inline std::atomic<std::int64_t> corrector_refine_hpn{0};        ///< HigherPrecisionNecessary from the refine corrector
inline std::atomic<std::int64_t> predictor_hpn{0};               ///< HigherPrecisionNecessary from the predictor
inline std::atomic<unsigned>     max_precision_seen{0};          ///< high-water working precision
inline std::atomic<unsigned>     max_digits_b{0};                ///< high-water DigitsB (the value that forces min_precision up)
inline std::atomic<unsigned>     max_startpoint_log10_affnorm{0}; ///< high-water log10||affine_solution|| in MHom start-point generation
inline std::atomic<unsigned>     max_startpoint_log10_Acond{0};   ///< high-water log10 cond(A) of the start-point linear system (trace-gated)

inline void reset()
{
	tracker_precision_increases = 0;
	endgame_refine_escalations  = 0;
	corrector_track_hpn         = 0;
	corrector_refine_hpn        = 0;
	predictor_hpn               = 0;
	max_precision_seen          = 0;
	max_digits_b                = 0;
	max_startpoint_log10_affnorm = 0;
	max_startpoint_log10_Acond   = 0;
}

inline void note_max(std::atomic<unsigned>& slot, unsigned v)
{
	unsigned cur = slot.load();
	while (v > cur && !slot.compare_exchange_weak(cur, v)) { /* retry */ }
}

inline void note_precision(unsigned p) { note_max(max_precision_seen, p); }
inline void note_digits_b(unsigned d)  { note_max(max_digits_b, d); }

// log10 magnitude helper (clamps non-finite/<=0 to 0); for the start-point conditioning probes.
inline void note_startpoint_affnorm(double affnorm)
{
	if (affnorm > 1.0) note_max(max_startpoint_log10_affnorm, static_cast<unsigned>(std::log10(affnorm)));
}
inline void note_startpoint_Acond(double cond)
{
	if (cond > 1.0) note_max(max_startpoint_log10_Acond, static_cast<unsigned>(std::log10(cond)));
}

} // namespace probe
} // namespace bertini
