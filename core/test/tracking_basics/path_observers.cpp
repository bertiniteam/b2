//This file is part of Bertini 2.
//
//path_observers.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//path_observers.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with path_observers.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire




#include <boost/test/unit_test.hpp>

#include "bertini2/trackers/amp_tracker.hpp"
#include "bertini2/trackers/observers.hpp"




extern double threshold_clearance_d;
extern bertini::mpfr_float threshold_clearance_mp;
extern unsigned TRACKING_TEST_MPFR_DEFAULT_DIGITS;




BOOST_AUTO_TEST_SUITE(AMP_tracker_basics)



using System = bertini::System;
using Variable = bertini::node::Variable;

using Var = std::shared_ptr<Variable>;

using VariableGroup = bertini::VariableGroup;


using dbl = std::complex<double>;
using mpfr = bertini::mpfr_complex;
using mpfr_float = bertini::mpfr_float;


template<typename NumT> using Vec = bertini::Vec<NumT>;
template<typename NumT> using Mat = bertini::Mat<NumT>;
using bertini::DefaultPrecision;


BOOST_AUTO_TEST_CASE(accumulate_single_path_square_root)
{
	DefaultPrecision(16);
	using namespace bertini::tracking;

	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var t = Variable::Make("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddFunction(x-t);
	sys.AddFunction(pow(y,2)-x);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;


	tracker.Setup(Predictor::Euler,
	              	1e-5,
					1e5,
					stepping_preferences,
					newton_preferences);

	tracker.PrecisionSetup(AMP);

	mpfr t_start(1);
	mpfr t_end(0);
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;


	AMPPathAccumulator<AMPTracker> path_accumulator;
	PrecisionAccumulator<AMPTracker> precision_accumulator;

	tracker.AddObserver(path_accumulator);
	tracker.AddObserver(precision_accumulator);

	start_point << mpfr(1), mpfr(1);
	[[maybe_unused]] bertini::SuccessCode tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);	
}




BOOST_AUTO_TEST_CASE(some_other_thing_square_root)
{
	DefaultPrecision(16);
	using namespace bertini::tracking;

	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var t = Variable::Make("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddFunction(x-t);
	sys.AddFunction(pow(y,2)-x);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;


	tracker.Setup(Predictor::Euler,
	              	1e-5,
					1e5,
					stepping_preferences,
					newton_preferences);

	tracker.PrecisionSetup(AMP);

	mpfr t_start(1);
	mpfr t_end(0);
	
	Vec<mpfr> start_point(2);
	Vec<mpfr> end_point;

	[[maybe_unused]] bertini::SuccessCode tracking_success;

	GoryDetailLogger<AMPTracker> tons_of_detail;

	tracker.AddObserver(tons_of_detail);

	start_point << mpfr(1), mpfr(1);
	tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

}



BOOST_AUTO_TEST_CASE(union_of_observers)
{
	DefaultPrecision(16);
	using namespace bertini::tracking;

	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var t = Variable::Make("t");

	System sys;

	VariableGroup v{x,y};

	sys.AddFunction(x-t);
	sys.AddFunction(pow(y,2)-x);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);

	auto AMP = bertini::tracking::AMPConfigFrom(sys);

	bertini::tracking::AMPTracker tracker(sys);


	SteppingConfig stepping_preferences;
	NewtonConfig newton_preferences;


	tracker.Setup(Predictor::Euler,
	              	1e-5,
					1e5,
					stepping_preferences,
					newton_preferences);

	tracker.PrecisionSetup(AMP);

	mpfr t_start(1);
	mpfr t_end(0);
	
	Vec<mpfr> start_point(2);
	start_point << mpfr(1), mpfr(1);

	Vec<mpfr> end_point;

	bertini::MultiObserver<AMPTracker, GoryDetailLogger> agglomeration;
	tracker.AddObserver(agglomeration);

	
	[[maybe_unused]] bertini::SuccessCode tracking_success = tracker.TrackPath(end_point,
	                  t_start, t_end, start_point);

}






// ---------------------------------------------------------------------------
// Observer-dispatch refactor (#255): self-unsubscribe via return value, and
// subscribing a new observer from inside Observe() (the "meta-observer" trick).
// ---------------------------------------------------------------------------

using bertini::AnyEvent;
using bertini::ObserveResult;
using bertini::tracking::AMPTracker;
using EmitterT = bertini::tracking::TrackerTraits<AMPTracker>::EventEmitterType;

// Counts every event it receives, and keeps observing forever.
template<class TrackerT>
struct CountingObserver : public bertini::Observer<TrackerT>
{
	int count = 0;
	ObserveResult Observe(AnyEvent const&) override
	{
		++count;
		return ObserveResult::KeepObserving;
	}
};

// Counts events, but asks to be unsubscribed after the very first one.
template<class TrackerT>
struct OneShotObserver : public bertini::Observer<TrackerT>
{
	int count = 0;
	ObserveResult Observe(AnyEvent const&) override
	{
		++count;
		return ObserveResult::Unsubscribe;
	}
};

// On TrackingStarted, attaches `child` to the emitting tracker; on TrackingEnded,
// detaches it.  Both happen from inside Observe(), exercising deferred mutation.
template<class TrackerT>
struct MetaObserver : public bertini::Observer<TrackerT>
{
	CountingObserver<TrackerT> child;
	bool saw_started = false;
	bool saw_ended   = false;

	ObserveResult Observe(AnyEvent const& e) override
	{
		if (auto p = dynamic_cast<const bertini::tracking::TrackingStarted<EmitterT>*>(&e))
		{
			saw_started = true;
			p->Get().AddObserver(child);   // deferred: child won't see THIS event
		}
		else if (auto p = dynamic_cast<const bertini::tracking::TrackingEnded<EmitterT>*>(&e))
		{
			saw_ended = true;
			p->Get().RemoveObserver(child);
		}
		return ObserveResult::KeepObserving;
	}
};


// builds the square-root system into `sys`, ready for a tracker.
static void BuildSquareRootSystem(System& sys)
{
	Var x = Variable::Make("x");
	Var y = Variable::Make("y");
	Var t = Variable::Make("t");

	VariableGroup v{x,y};
	sys.AddFunction(x-t);
	sys.AddFunction(pow(y,2)-x);
	sys.AddPathVariable(t);
	sys.AddVariableGroup(v);
}


BOOST_AUTO_TEST_CASE(self_unsubscribe_via_return_value)
{
	DefaultPrecision(16);
	using namespace bertini::tracking;

	System sys;
	BuildSquareRootSystem(sys);
	AMPTracker tracker(sys);
	tracker.Setup(Predictor::Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig());
	tracker.PrecisionSetup(AMPConfigFrom(sys));

	OneShotObserver<AMPTracker> one_shot;
	tracker.AddObserver(one_shot);

	Vec<mpfr> start_point(2);
	start_point << mpfr(1), mpfr(1);
	Vec<mpfr> end_point;
	[[maybe_unused]] auto code = tracker.TrackPath(end_point, mpfr(1), mpfr(0), start_point);

	// it returned Unsubscribe on the first event, so it must have been dropped
	// before the second event was ever emitted.
	BOOST_CHECK_EQUAL(one_shot.count, 1);
}


BOOST_AUTO_TEST_CASE(meta_observer_attaches_child_mid_dispatch)
{
	DefaultPrecision(16);
	using namespace bertini::tracking;

	System sys;
	BuildSquareRootSystem(sys);
	AMPTracker tracker(sys);
	tracker.Setup(Predictor::Euler, 1e-5, 1e5, SteppingConfig(), NewtonConfig());
	tracker.PrecisionSetup(AMPConfigFrom(sys));

	MetaObserver<AMPTracker> meta;
	tracker.AddObserver(meta);

	Vec<mpfr> start_point(2);
	start_point << mpfr(1), mpfr(1);
	Vec<mpfr> end_point;
	[[maybe_unused]] auto code = tracker.TrackPath(end_point, mpfr(1), mpfr(0), start_point);

	BOOST_CHECK(meta.saw_started);
	BOOST_CHECK(meta.saw_ended);
	// the child was attached during TrackingStarted and saw the events that
	// followed (steps, precision changes, ...) up to and including TrackingEnded.
	BOOST_CHECK_GT(meta.child.count, 0);
}


BOOST_AUTO_TEST_SUITE_END()



