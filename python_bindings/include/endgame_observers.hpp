//This file is part of Bertini 2.
//
//python/endgame_observers.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/endgame_observers.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/endgame_observers.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  silviana amethyst
//  UWEC
//  Fall 2017, Spring 2018
//
//
//  python/endgame_observers.hpp:  source file for exposing endgame observers to python.

#pragma once
#include "python_common.hpp"
#include "endgame_export.hpp"

#include <bertini2/endgames/observers.hpp>


namespace bertini{
	namespace python{


void ExportEndgameObservers();


using namespace bertini::endgame;


template<typename ObsT>
struct EndgameObserverVisitor: public def_visitor<EndgameObserverVisitor<ObsT> >
{
	friend class ::boost::python::def_visitor_access;

public:

	template<class PyClass>
	void visit(PyClass& /*cl*/) const{
	}
};


/**
\brief Copy a container of space points into a Python list.

Boost.Python has no converter for std::vector<Vec<T>>, but eigenpy converts each
Vec<T> on its own, so the bucket is handed over one element at a time.
*/
template<typename ContT>
boost::python::list VectorOfPointsToList(ContT const& c)
{
	boost::python::list out;
	for (auto const& v : c)
		out.append(v);
	return out;
}

/**
\brief Copy a container of scalars into a Python list.
*/
template<typename ContT>
boost::python::list VectorOfScalarsToList(ContT const& c)
{
	boost::python::list out;
	for (auto const& v : c)
		out.append(v);
	return out;
}


/**
\brief Exposes the three buckets of a SampleSequenceCollector as Python lists.

Each bucket is copied out on access -- the collector owns its vectors and keeps
collecting, so handing Python a view would alias a container that grows.
*/
template<typename ObsT>
struct SampleSequenceVisitor: public def_visitor<SampleSequenceVisitor<ObsT> >
{
	friend class ::boost::python::def_visitor_access;

public:

	template<class PyClass>
	void visit(PyClass& cl) const {
		cl
		.def("path_samples", &PathSamples, arg("self"),
			"The sequence: points ON the path, at geometrically shrinking times, approaching the root.  "
			"These are what a trend read wants -- a quantity that vanishes at the root traces a power "
			"law against the matching path_times, and one that does not is flat.")
		.def("path_times", &PathTimes, arg("self"),
			"The time at which each path sample was computed.  Same length as path_samples.")
		.def("circle_samples", &CircleSamples, arg("self"),
			"Circle-track points (Cauchy only).  These sit at CONSTANT |t| and do not approach the "
			"root, so they are NOT part of the sequence -- their mean is what becomes an approximation.  "
			"Kept separate for loop diagnostics.")
		.def("circle_times", &CircleTimes, arg("self"),
			"The time of each circle sample.  Same length as circle_samples.")
		.def("approximations", &Approximations, arg("self"),
			"Successive approximations OF the root, as the endgame produced them -- estimates of the "
			"root rather than points on the path.")
		.def("approximation_times", &ApproximationTimes, arg("self"),
			"The endgame's latest time at each approximation.")
		.def("approximation_errors", &ApproximationErrors, arg("self"),
			"The error the endgame reported at each approximation -- the quantity its own convergence "
			"test compares against the final tolerance.")
		.def("cycle_numbers", &CycleNumbers, arg("self"),
			"The cycle number the endgame held at each approximation.  The exponent in the power law: "
			"distance to the root falls as |t|^(1/c).")
		.def("advance_times", &AdvanceTimes, arg("self"),
			"Times at which TimeAdvanced fired.  That event carries no payload of its own, and is "
			"emitted by the Cauchy endgame only.")
		.def("num_samples", &ObsT::NumSamples, arg("self"),
			"The number of samples collected -- the length of the sequence.")
		.def("num_precision_increases", &NumPrecisionIncreases, arg("self"),
			"How many times the endgame raised its working precision, over every run observed.  When "
			"that happens before the first approximation exists, the endgame recomputes its sample "
			"window at the new precision and the superseded samples are dropped from the sequence, so "
			"path_times keep marching toward the target within each run.")
		.def("num_runs", &ObsT::NumRuns, arg("self"),
			"The number of endgame runs observed -- the number of paths, when attached to a solver's "
			"endgame via get_endgame().")
		.def("run_path_starts", &RunPathStarts, arg("self"),
			"Index into path_samples where each run's samples begin, so one path's sequence can be told "
			"from the next: run j owns path_samples[run_path_starts[j] : run_path_starts[j+1]].")
		.def("run_circle_starts", &RunCircleStarts, arg("self"),
			"Index into circle_samples where each run's circle points begin.")
		.def("run_approx_starts", &RunApproxStarts, arg("self"),
			"Index into approximations where each run's approximations begin.")
		.def("clear", &ObsT::Clear, arg("self"),
			"Forget everything collected so far, so one collector can be reused across paths.")
		;
	}

private:
	static size_t NumPrecisionIncreases(ObsT const& self)            { return self.num_precision_increases; }
	static boost::python::list PathSamples(ObsT const& self)         { return VectorOfPointsToList(self.path_samples); }
	static boost::python::list PathTimes(ObsT const& self)           { return VectorOfScalarsToList(self.path_times); }
	static boost::python::list CircleSamples(ObsT const& self)       { return VectorOfPointsToList(self.circle_samples); }
	static boost::python::list CircleTimes(ObsT const& self)         { return VectorOfScalarsToList(self.circle_times); }
	static boost::python::list Approximations(ObsT const& self)      { return VectorOfPointsToList(self.approximations); }
	static boost::python::list ApproximationTimes(ObsT const& self)  { return VectorOfScalarsToList(self.approximation_times); }
	static boost::python::list ApproximationErrors(ObsT const& self) { return VectorOfScalarsToList(self.approximation_errors); }
	static boost::python::list CycleNumbers(ObsT const& self)        { return VectorOfScalarsToList(self.cycle_numbers); }
	static boost::python::list AdvanceTimes(ObsT const& self)        { return VectorOfScalarsToList(self.advance_times); }
	static boost::python::list RunPathStarts(ObsT const& self)       { return VectorOfScalarsToList(self.run_path_starts); }
	static boost::python::list RunCircleStarts(ObsT const& self)     { return VectorOfScalarsToList(self.run_circle_starts); }
	static boost::python::list RunApproxStarts(ObsT const& self)     { return VectorOfScalarsToList(self.run_approx_starts); }
};



}} // namespaces
