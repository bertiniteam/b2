//This file is part of Bertini 2.
//
//endgames/observers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//endgames/observers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with endgames/observers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file endgames/observers.hpp

\brief Contains the endgames/observers base types
*/

#pragma once

#include "bertini2/endgames/events.hpp"
#include "bertini2/logging.hpp"
#include "bertini2/detail/observer.hpp"

#include <boost/type_index.hpp>

#include <vector>

namespace bertini {

    namespace endgame{


/**
\brief Logs the endgame run, with gory detail.

\ingroup observer
*/
template <typename EndgameT>
struct GoryDetailLogger : public Observer<EndgameT>
{BOOST_TYPE_INDEX_REGISTER_CLASS

using EmitterT = EndgameT;  ///< The event-emitter type.
using BCT = typename EndgameT::BaseComplexT;  ///< The complex number type.

virtual ~GoryDetailLogger() = default;

virtual ObserveResult Observe(AnyEvent const& e) override
{
    if(auto p = dynamic_cast<const TimeAdvanced<EmitterT>*>(&e))
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "time advanced " << p->Get().LatestTime();
    }

    else if (auto p = dynamic_cast<const SampleRefined<EmitterT>*>(&e))
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "refined a sample, huzzah";
    }

    else if (auto p = dynamic_cast<const CircleAdvanced<EmitterT>*>(&e))
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "advanced around the circle, to " << p->NewSample()<< " at time " << p->NewTime();
    }

    else if (auto p = dynamic_cast<const ClosedLoop<EmitterT>*>(&e))
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "closed a loop, cycle number " << p->Get().CycleNumber();
    }
    else if (auto p = dynamic_cast<const ApproximatedRoot<EmitterT>*>(&e))
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "approximated the target root.  approximation " << p->Get().template FinalApproximation<BCT>() << " with error " << p->Get().ApproximateError();
    }

    else if (auto p = dynamic_cast<const PrecisionChanged<AMPEndgame>*>(&e))
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "precision changed from  " << p->Previous() << " to " << p->Next();
    }

    else if (auto p = dynamic_cast<const InEGOperatingZone<EmitterT>*>(&e))
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "made it to the endgame operating zone at time " << p->Get().LatestTime();
    }

    else if(auto p = dynamic_cast<const Converged<EmitterT>*>(&e))
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "converged at time " << p->Get().LatestTime() << " with result " << p->Get().template FinalApproximation<BCT>() << " and residual " << p->Get().ApproximateError();
    }
    else if (auto p = dynamic_cast<const Initializing<EmitterT>*>(&e))
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "starting running " << boost::typeindex::type_id<EmitterT>().pretty_name();
    }
    else
    {
        BOOST_LOG_TRIVIAL(severity_level::debug) << "unprogrammed response for event of type " << boost::typeindex::type_id_runtime(e).pretty_name();
    }

    return ObserveResult::KeepObserving;
}

}; // gory detail



/**
\brief Counts endgame events and captures their payloads, for tests and diagnostics.

Records how many of each event type were delivered, and captures the most-recent CircleAdvanced
point/time and the Converged approximation -- reading each payload (and the state accessors LatestTime /
FinalApproximation / ApproximateError) the way a real observer would.  This deliberately exercises the
full event-delivery path, including the numeric-type conversion the adaptive-numeric-type endgame
performs when it emits from its hardware-complex_dbl fast lane, and the slot-aware state accessors.

\ingroup observer
*/
template <typename EndgameT>
struct EventRecorder : public Observer<EndgameT>
{BOOST_TYPE_INDEX_REGISTER_CLASS

using EmitterT = EndgameT;                     ///< The endgame type emitting the observed events.
using BCT = typename EndgameT::BaseComplexT;    ///< The boundary complex type of the observed endgame.

unsigned num_events            = 0;            ///< Total number of events observed.
unsigned num_time_advanced     = 0;            ///< Number of TimeAdvanced events observed.
unsigned num_sample_refined    = 0;            ///< Number of SampleRefined events observed.
unsigned num_circle_advanced   = 0;            ///< Number of circle-advanced events observed.
unsigned num_closed_loop       = 0;            ///< Number of closed-loop events observed.
unsigned num_approximated_root = 0;            ///< Number of approximated-root events observed.
unsigned num_in_eg_zone        = 0;            ///< Number of in-endgame-zone events observed.
unsigned num_converged         = 0;            ///< Number of converged events observed.
unsigned num_precision_changed = 0;            ///< Number of precision-changed events observed.

Vec<BCT> last_circle_point;                    ///< The most-recent circle sample point captured from an event.
BCT      last_circle_time;                     ///< The most-recent circle sample time captured from an event.
Vec<BCT> converged_point;                      ///< The converged root point captured from a converged event.

virtual ObserveResult Observe(AnyEvent const& e) override
{
    ++num_events;
    if (auto p = dynamic_cast<const TimeAdvanced<EmitterT>*>(&e))
    { ++num_time_advanced; (void)p->Get().LatestTime(); }

    else if (dynamic_cast<const SampleRefined<EmitterT>*>(&e))
    { ++num_sample_refined; }

    else if (auto p = dynamic_cast<const CircleAdvanced<EmitterT>*>(&e))
    { ++num_circle_advanced; last_circle_point = p->NewSample(); last_circle_time = p->NewTime(); }

    else if (dynamic_cast<const ClosedLoop<EmitterT>*>(&e))
    { ++num_closed_loop; }

    else if (auto p = dynamic_cast<const ApproximatedRoot<EmitterT>*>(&e))
    { ++num_approximated_root; (void)p->Get().template FinalApproximation<BCT>(); (void)p->Get().ApproximateError(); }

    else if (dynamic_cast<const InEGOperatingZone<EmitterT>*>(&e))
    { ++num_in_eg_zone; }

    else if (auto p = dynamic_cast<const Converged<EmitterT>*>(&e))
    { ++num_converged; converged_point = p->Get().template FinalApproximation<BCT>(); (void)p->Get().LatestTime(); }

    else if (dynamic_cast<const PrecisionChanged<AMPEndgame>*>(&e))
    { ++num_precision_changed; }

    return ObserveResult::KeepObserving;
}

}; // EventRecorder



/**
\brief Collects an endgame's approach to the root as a time-indexed sequence, for callers
that need to watch a quantity BEHAVE as the point improves rather than judge it at a
single point.

A single spectrum cannot separate a genuinely tiny singular value from a perturbation
artifact -- the two are numerically identical at one point.  What separates them is the
trend: a truly-zero value is bounded by the point's own error, so it tracks the distance
to the root down, while a genuinely nonzero one plateaus.  Reading that trend needs the
approach itself, not just its final answer, which is what this collects.

THREE BUCKETS, kept apart because they are different evidence:

- ``path_samples`` -- the samples, from ComputedSamplePoint.  Points ON the path at
  geometrically shrinking times, approaching the root.  These are the sequence: distance to
  the root falls as |t|^(1/c) with c the cycle number, so a quantity that vanishes at the
  root traces a power law against these times, and one that does not is flat.
- ``circle_samples`` -- from CircleAdvanced, Cauchy only.  These sit at CONSTANT |t| and
  do not approach the root, so they are not part of the sequence; their mean is what becomes an
  approximation.  Kept for loop diagnostics, never mixed into the sequence.
- ``approximations`` -- from ApproximatedRoot, with the error and cycle number the
  endgame reported at each.  Estimates OF the root rather than points on the path.

``advance_times`` records TimeAdvanced, which carries no payload of its own -- it marks
that the approach axis moved, and is emitted by the Cauchy endgame only.

Attach one of these per path.  It is an ordinary C++ observer, so the observable filters
by event type before any virtual call, and a run that emits thousands of events pays only
for the handful this asks for.

\ingroup observer
*/
template <typename EndgameT>
struct SampleSequenceCollector : public Observer<EndgameT>
{BOOST_TYPE_INDEX_REGISTER_CLASS

using EmitterT = EndgameT;                    ///< The endgame type emitting the observed events.
using BCT = typename EndgameT::BaseComplexT;  ///< The boundary complex type of the observed endgame.

std::vector<Vec<BCT>> path_samples;      ///< Points on the path, from ComputedSamplePoint -- the sequence.
std::vector<BCT>      path_times;        ///< The time at which each path sample was computed.

std::vector<Vec<BCT>> circle_samples;    ///< Circle-track points, from CircleAdvanced (Cauchy only).
std::vector<BCT>      circle_times;      ///< The time of each circle sample.

std::vector<Vec<BCT>> approximations;      ///< Root approximations, from ApproximatedRoot.
std::vector<BCT>      approximation_times; ///< The endgame's latest time at each approximation.
std::vector<NumErrorT> approximation_errors; ///< The endgame's reported error at each approximation.
std::vector<unsigned>  cycle_numbers;      ///< The endgame's cycle number at each approximation.

std::vector<BCT> advance_times;          ///< Times at which TimeAdvanced fired (that event carries no payload).

// Run boundaries.  One entry per endgame Run, holding the size each bucket had when that
// run began -- so a collector attached to a whole solve can still tell one path from the
// next.  Samples of run j are path_samples[run_path_starts[j] .. run_path_starts[j+1]).
std::vector<size_t> run_path_starts;     ///< Index into path_samples where each run's samples begin.
std::vector<size_t> run_circle_starts;   ///< Index into circle_samples where each run's circle points begin.
std::vector<size_t> run_approx_starts;   ///< Index into approximations where each run's approximations begin.

size_t num_precision_increases = 0;   ///< How many times the endgame raised its working precision (PrecisionChanged), over every run observed.  When that happens before the first approximation the sample window is recomputed at the new precision and the superseded samples are dropped from the sequence (SamplesRecomputedAtHigherPrecision).

/// \brief Forget everything collected so far, so one collector can be reused across paths.
void Clear()
{
    path_samples.clear();          path_times.clear();
    circle_samples.clear();        circle_times.clear();
    approximations.clear();        approximation_times.clear();
    approximation_errors.clear();  cycle_numbers.clear();
    advance_times.clear();
    run_path_starts.clear();       run_circle_starts.clear();
    run_approx_starts.clear();
    num_precision_increases = 0;
}

/// \return The number of endgame runs observed -- the number of paths, when attached to a solver.
size_t NumRuns() const { return run_path_starts.size(); }

/// \return The number of samples collected -- the length of the sequence.
size_t NumSamples() const { return path_samples.size(); }

virtual ObserveResult Observe(AnyEvent const& e) override
{
    if (auto p = dynamic_cast<const ComputedSamplePoint<EmitterT>*>(&e))
    {
        path_samples.push_back(p->NewSample());
        path_times.push_back(p->NewTime());
    }

    else if (auto p = dynamic_cast<const CircleAdvanced<EmitterT>*>(&e))
    {
        circle_samples.push_back(p->NewSample());
        circle_times.push_back(p->NewTime());
    }

    else if (auto p = dynamic_cast<const ApproximatedRoot<EmitterT>*>(&e))
    {
        approximations.push_back(p->Get().template FinalApproximation<BCT>());
        approximation_times.push_back(p->Get().LatestTime());
        approximation_errors.push_back(p->Get().ApproximateError());
        cycle_numbers.push_back(p->Get().CycleNumber());
    }

    else if (auto p = dynamic_cast<const TimeAdvanced<EmitterT>*>(&e))
    {
        advance_times.push_back(p->Get().LatestTime());
    }

    else if (dynamic_cast<const Initializing<EmitterT>*>(&e))
    {
        run_path_starts.push_back(path_samples.size());
        run_circle_starts.push_back(circle_samples.size());
        run_approx_starts.push_back(approximations.size());
    }

    else if (dynamic_cast<const PrecisionChanged<EmitterT>*>(&e))
    {
        ++num_precision_increases;
    }

    else if (dynamic_cast<const SamplesRecomputedAtHigherPrecision<EmitterT>*>(&e))
    {
        // The endgame needed a higher precision before it had its first approximation, and
        // tracks its sample window again at the new precision.  The samples it announced at
        // the lower precision are superseded: drop them as the endgame does, so the sequence
        // holds one approach at one precision and its times keep marching toward the target.
        if (!run_path_starts.empty())
        {
            path_samples.resize(run_path_starts.back());
            path_times.resize(run_path_starts.back());
            circle_samples.resize(run_circle_starts.back());
            circle_times.resize(run_circle_starts.back());
            approximations.resize(run_approx_starts.back());
            approximation_times.resize(run_approx_starts.back());
            approximation_errors.resize(run_approx_starts.back());
            cycle_numbers.resize(run_approx_starts.back());
        }
    }

    return ObserveResult::KeepObserving;
}

}; // SampleSequenceCollector


    } //re: namespace endgames
}// re: namespace bertini
