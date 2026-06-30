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



/**
\brief Runtime options selecting which ZeroDim algorithm to instantiate.

Holds the runtime (enum) choices -- start system, tracker, and endgame -- that the
ZeroDimSpecify* switch chain resolves into concrete compile-time template parameters.
*/
struct ZeroDimRT
{
	// INTERIM default = RootsOfUnity (see policies.hpp): the linear-product TotalDegree is the
	// eventual default but currently stalls the Cauchy endgame on harder systems, so the safe
	// default stays roots of unity until that is fixed.
	type::Start start = type::Start::RootsOfUnity;    ///< Which start system to use.
	type::Tracker tracker = type::Tracker::Adaptive;  ///< Which path tracker to use.
	type::Endgame endgame = type::Endgame::Cauchy;    ///< Which endgame to use.
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
// start_system::MakeStartFactory<StartType>() (see start_base.hpp); ZeroDim downstream is a
// single type that holds the start system polymorphically.  Add a start system => add a
// case in ZeroDimSpecifyStart, no new ZeroDim instantiation.

/**
\brief Instantiate the fully-specified ZeroDim solver -- the end of the switch chain.

By this point every algorithm choice is a compile-time template parameter, so the concrete
ZeroDimSolver can be constructed.  The blackbox always clones and builds its own system.

\tparam TrackerType The resolved path-tracker type.
\tparam EndgameType The resolved endgame type.
\tparam ConstTs The forwarded construction-argument types.
\param ts The construction arguments forwarded to the solver (target system and start-system factory).
\return An owning handle to the constructed solver, type-erased as AnyZeroDim.
*/
template <typename TrackerType, typename EndgameType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyComplete(ConstTs const& ...ts)
{
	// the blackbox always clones + builds (ZeroDimSolver); ts... = (target, start_factory).
	return std::make_unique<
			algorithm::ZeroDimSolver<TrackerType, EndgameType, System>
			>(ts...);
}

/**
\brief Resolve the runtime endgame choice (rt.endgame) into a compile-time endgame type, then continue.

\tparam TrackerType The already-resolved path-tracker type.
\tparam ConstTs The forwarded construction-argument types.
\param rt The runtime options carrying the endgame selection.
\param ts The construction arguments forwarded down the chain.
\return An owning handle to the constructed solver, type-erased as AnyZeroDim.
*/
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

/**
\brief Resolve the runtime tracker choice (rt.tracker) into a compile-time tracker type, then continue.

\tparam ConstTs The forwarded construction-argument types.
\param rt The runtime options carrying the tracker selection.
\param ts The construction arguments forwarded down the chain.
\return An owning handle to the constructed solver, type-erased as AnyZeroDim.
*/
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

/**
\brief Resolve the runtime start-system choice (rt.start) into a start-system factory, then continue.

Appends the chosen start-system factory to the argument pack before continuing down the chain
(ZeroDimSolver consumes the target system and that factory).

\tparam ConstTs The forwarded construction-argument types.
\param rt The runtime options carrying the start-system selection.
\param ts The construction arguments forwarded down the chain.
\return An owning handle to the constructed solver, type-erased as AnyZeroDim.
*/
template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyStart(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	// append the start-system factory to the argument pack; ZeroDimSolver consumes (target, factory).
	switch (rt.start)
	{
		case type::Start::TotalDegree:
			return ZeroDimSpecifyTracker(rt, ts..., start_system::MakeStartFactory<start_system::TotalDegree>());
		case type::Start::RootsOfUnity:
			return ZeroDimSpecifyTracker(rt, ts..., start_system::MakeStartFactory<start_system::RootsOfUnity>());
		case type::Start::MHom:
			return ZeroDimSpecifyTracker(rt, ts..., start_system::MakeStartFactory<start_system::MHomogeneous>());
		case type::Start::User:
			throw std::runtime_error("trying to use generic zero dim with user homotopy.  use the specific UserBlaBla instead");
	}
	throw std::runtime_error("unrecognized start system type in ZeroDimSpecifyStart");
}

/**
\brief Build a ZeroDim solver from runtime options -- the entry point of the switch chain.

Enters the ZeroDimSpecify* chain, which resolves each runtime option (start, tracker, endgame)
into the corresponding compile-time template parameter and constructs the solver.

\tparam ConstTs The forwarded construction-argument types.
\param rt The runtime options selecting start system, tracker, and endgame.
\param ts The construction arguments (the target system, ...).
\return An owning handle to the constructed solver, type-erased as AnyZeroDim.
*/
template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> MakeZeroDim(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	return ZeroDimSpecifyStart(rt, ts...);
}


} //ns blackbox
} //ns bertini
