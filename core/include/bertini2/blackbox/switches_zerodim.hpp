//This file is part of Bertini 2.
//
//bertini2/blackbox/switches_zerodim.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/blackbox/switches_zerodim.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/blackbox/switches_zerodim.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.
//
// silviana amethyst, university of wisconsin-eau claire

/**
\file bertini2/blackbox/switches_zerodim.hpp 

\brief A sequence of switches for getting a particular instantiation of a ZeroDim algorithm, based on runtime options.
*/



#pragma once


#include "bertini2/system.hpp"
#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/endgames.hpp"

#include "bertini2/blackbox/config.hpp"


namespace bertini{
namespace blackbox{



struct ZeroDimRT
{
	// INTERIM default = RootsOfUnity (see policies.hpp): the linear-product TotalDegree is the
	// eventual default but currently stalls the Cauchy endgame on harder systems, so the safe
	// default stays roots of unity until that is fixed.
	type::Start start = type::Start::RootsOfUnity;
	type::Tracker tracker = type::Tracker::Adaptive;
	type::Endgame endgame = type::Endgame::Cauchy;
};


/**
\brief Infer which start system to use from the target system's variable-group structure.

This replicates classic Bertini, which chooses the start system from how the user
groups the variables rather than from a dedicated setting:

- a single affine variable group, with no homogeneous variable groups -> total degree
  (the 1-homogeneous Bezout start system);
- anything else with grouping -- two or more variable groups, or one or more
  homogeneous variable groups -- -> multihomogeneous, using that partition.

User-defined homotopies are not inferred here; they come with their own start system.
*/
inline type::Start InferStartType(System const& sys)
{
	// INTERIM: a single affine group infers RootsOfUnity (the safe default).  The eventual choice
	// here is TotalDegree (general position), gated on the Cauchy-endgame fix; see policies.hpp.
	if (sys.NumVariableGroups() == 1 && sys.NumHomVariableGroups() == 0)
		return type::Start::RootsOfUnity;
	else
		return type::Start::MHom;
}


// The concrete start-system type appears only at the construction site, via
// policy::MakeStartFactory<StartType>() (see common/policies.hpp); ZeroDim downstream is a
// single type that holds the start system polymorphically.  Add a start system => add a
// case in ZeroDimSpecifyStart, no new ZeroDim instantiation.

template <typename TrackerType, typename EndgameType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyComplete(ConstTs const& ...ts)
{
	// the blackbox always clones (default CloneGiven policy); ts... = (target, start_factory).
	return std::make_unique<
			algorithm::ZeroDim<TrackerType, EndgameType, System>
			>(ts...);
}

template <typename TrackerType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyEndgame(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.endgame)
	{
		case type::Endgame::PowerSeries:
			// honor the requested endgame!  until 2026-06-12 this hardcoded Cauchy.
			return ZeroDimSpecifyComplete<TrackerType,
					typename endgame::EndgameSelector<TrackerType>::PSEG>(ts...);

		case type::Endgame::Cauchy:
			return ZeroDimSpecifyComplete<TrackerType,
					typename endgame::EndgameSelector<TrackerType>::Cauchy>(ts...);
	}
	throw std::runtime_error("unrecognized endgame type in ZeroDimSpecifyEndgame");
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyTracker(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.tracker)
	{
		case type::Tracker::FixedDouble:
			return ZeroDimSpecifyEndgame<tracking::DoublePrecisionTracker>(rt, ts...);
		case type::Tracker::FixedMultiple:
			return ZeroDimSpecifyEndgame<tracking::MultiplePrecisionTracker>(rt, ts...);
		case type::Tracker::Adaptive:
			return ZeroDimSpecifyEndgame<tracking::AMPTracker>(rt, ts...);
	}
	throw std::runtime_error("unrecognized tracker type in ZeroDimSpecifyTracker");
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyStart(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	// append the start-system factory to the argument pack; CloneGiven consumes (target, factory).
	switch (rt.start)
	{
		case type::Start::TotalDegree:
			return ZeroDimSpecifyTracker(rt, ts..., policy::MakeStartFactory<start_system::TotalDegree>());
		case type::Start::RootsOfUnity:
			return ZeroDimSpecifyTracker(rt, ts..., policy::MakeStartFactory<start_system::RootsOfUnity>());
		case type::Start::MHom:
			return ZeroDimSpecifyTracker(rt, ts..., policy::MakeStartFactory<start_system::MHomogeneous>());
		case type::Start::User:
			throw std::runtime_error("trying to use generic zero dim with user homotopy.  use the specific UserBlaBla instead");
	}
	throw std::runtime_error("unrecognized start system type in ZeroDimSpecifyStart");
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> MakeZeroDim(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	return ZeroDimSpecifyStart(rt, ts...);
}


} //ns blackbox
} //ns bertini
