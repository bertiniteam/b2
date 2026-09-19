//This file is part of Bertini 2.
//
//operating_zone_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//operating_zone_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with operating_zone_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file InEGOperatingZone means the same thing in both endgames (b2#402).

The event announces that the path has reached the asymptotic regime, where the Puiseux model
the endgame is built on dominates.  That is what makes the samples usable quantitatively -- a
consumer fitting `sigma ~ |t|^(1/c)` needs to know which rungs lie on the power law.

Cauchy emitted it once, after its c/k stabilization loop.  The power series endgame emitted it
at the end of every successful advance, with no test at all, so there it meant "advanced" and a
consumer trusting it got a silently wrong answer under one flavor and a right one under the
other.  Both now run the same test, which lives on the base endgame because it reads only the
geometric approach samples that every flavor keeps.
*/

#include <boost/test/unit_test.hpp>

#include "bertini2/system/start_systems.hpp"
#include "bertini2/num_traits.hpp"
#include "bertini2/endgames/fixed_prec_endgame.hpp"
#include "bertini2/endgames/powerseries.hpp"
#include "bertini2/endgames/cauchy.hpp"
#include "bertini2/endgames/observers.hpp"

using namespace bertini;
using Variable = node::Variable;
using Var = std::shared_ptr<Variable>;

namespace {

using TrackerType = tracking::DoublePrecisionTracker;

/// Counts how many times an endgame says the path entered the operating zone.
template<typename EGT>
struct ZoneCounter : public Observer<EGT>
{
    using Emitter = EGT;
    unsigned entries = 0;  ///< InEGOperatingZone events seen.
    unsigned samples = 0;  ///< Sample points seen, for scale: the zone event must not track this.

    ObserveResult Observe(AnyEvent const& e) override
    {
        if (dynamic_cast<const endgame::InEGOperatingZone<Emitter>*>(&e))
            ++entries;
        if (dynamic_cast<const endgame::ComputedSamplePoint<Emitter>*>(&e))
            ++samples;
        return ObserveResult::KeepObserving;
    }
};

// x = 1 is a triple root at t = 0; a cycle number worth estimating, so the c/k estimates have
// something to settle to.
System TripleRoot()
{
    System sys;
    Var x = Variable::Make("x"), t = Variable::Make("t");
    sys.AddFunction( pow(x-1,3)*(1-t) + (pow(x,3)+1)*t );
    sys.AddVariableGroup(VariableGroup{x});
    sys.AddPathVariable(t);
    return sys;
}

template<typename EGT>
std::pair<unsigned, unsigned> RunAndCount()
{
    DefaultPrecision(16);
    auto sys = TripleRoot();

    TrackerType tracker(sys);
    tracking::SteppingConfig stepping;
    tracking::NewtonConfig newton;
    tracker.Setup(tracking::Predictor::HeunEuler, 1e-6, 1e5, stepping, newton);
    tracker.PrecisionSetup(tracking::FixedPrecisionConfig(sys));

    complex_dbl boundary_time = complex_dbl(0.1);
    bertini::Vec<complex_dbl> boundary_point(1);
    boundary_point << complex_dbl(5.000000000000001e-01, 9.084258952712920e-17);

    EGT endgame(tracker, endgame::EndgameConfig(), endgame::SecurityConfig());
    endgame.SetBoundaryTime(boundary_time);

    ZoneCounter<EGT> counter;
    endgame.AddObserver(counter);
    endgame.Run(boundary_point);
    endgame.RemoveObserver(counter);

    return {counter.entries, counter.samples};
}

} // unnamed namespace


BOOST_AUTO_TEST_SUITE(endgame_operating_zone)

// Entering the zone happens once, so the event fires once.  Under the old power-series code this
// counted one per advance -- the assertion that fails loudest if the criterion is removed.
BOOST_AUTO_TEST_CASE(the_power_series_endgame_announces_the_zone_once)
{
    auto const [entries, samples] = RunAndCount<endgame::EndgameSelector<TrackerType>::PSEG>();
    BOOST_REQUIRE(samples > 1);            // it really did advance more than once
    BOOST_CHECK_EQUAL(entries, 1u);
}


BOOST_AUTO_TEST_CASE(the_cauchy_endgame_announces_the_zone_once)
{
    auto const [entries, samples] = RunAndCount<endgame::EndgameSelector<TrackerType>::Cauchy>();
    (void)samples;
    BOOST_CHECK_EQUAL(entries, 1u);
}


// The criterion, on its own: estimates that agree are in the zone, estimates that do not are not.
// Both flavors reach this through the base, so neither can drift to a private notion of the zone.
BOOST_AUTO_TEST_CASE(the_criterion_is_the_settling_of_the_c_over_k_estimates)
{
    DefaultPrecision(16);
    auto sys = TripleRoot();
    TrackerType tracker(sys);
    tracking::SteppingConfig stepping;
    tracking::NewtonConfig newton;
    tracker.Setup(tracking::Predictor::HeunEuler, 1e-6, 1e5, stepping, newton);
    tracker.PrecisionSetup(tracking::FixedPrecisionConfig(sys));

    endgame::EndgameConfig eg_settings;
    using EGT = endgame::EndgameSelector<TrackerType>::PSEG;
    EGT endgame(tracker, eg_settings, endgame::SecurityConfig());

    // the default asks for three estimates agreeing to within 3/4
    std::deque<double> settled{2.0, 2.0, 2.0, 2.0};
    BOOST_CHECK(endgame.CheckForCOverKStabilization(settled));

    std::deque<double> wandering{1.0, 4.0, 1.0, 4.0};
    BOOST_CHECK(!endgame.CheckForCOverKStabilization(wandering));

    // and it is the user's threshold, not a hard-coded one
    std::deque<double> loose{1.0, 1.5, 1.0, 1.5};
    BOOST_CHECK(!endgame.CheckForCOverKStabilization(loose));
    eg_settings.minimum_for_c_over_k_stabilization = 0.5;
    endgame.Set(eg_settings);
    BOOST_CHECK(endgame.CheckForCOverKStabilization(loose));
}

BOOST_AUTO_TEST_SUITE_END()
