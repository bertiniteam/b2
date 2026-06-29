//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/numerical_irreducible_decomposition.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/numerical_irreducible_decomposition.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/numerical_irreducible_decomposition.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/nag_algorithms/numerical_irreducible_decomposition.hpp 

\brief Provides the NID algorithms for Bertini2.
*/



#pragma once

#include "bertini2/num_traits.hpp"

#include "bertini2/detail/visitable.hpp"
#include "bertini2/tracking.hpp"

#include "bertini2/detail/configured.hpp"
#include "bertini2/detail/observable.hpp"

#include "bertini2/nag_algorithms/common/algorithm_base.hpp"
#include "bertini2/nag_algorithms/common/config.hpp"

#include "bertini2/nag_datatypes/numerical_irreducible_decomposition.hpp"

#include <stdexcept>


namespace bertini {

	namespace algorithm {


/**
forward declare of the NumericalIrreducibleDecomposition algorithm.

Unlike ZeroDim, NID has no start system -- the regenerative cascade is a
fundamentally different way of solving a polynomial system -- so it carries no
StartSystem template parameter.  It owns a cloned target system directly (no
system-management policy: this is placeholder scaffolding whose Solve() throws).
*/
template<	typename TrackerType, typename EndgameType,
			typename SystemType = System >
struct NumericalIrreducibleDecomposition;



/**
specify the traits for the algorithm.  this is why we need the forward declare.

The NeededConfigs typelist is what drives the (reusable) Python config interface:
the ConfiguredVisitor reflects over exactly these types.
*/
template<	typename TrackerType, typename EndgameType,
			typename SystemType >
struct AlgoTraits< NumericalIrreducibleDecomposition<TrackerType, EndgameType, SystemType> >
{
	using BaseRealT    = typename tracking::TrackerTraits<TrackerType>::BaseRealT;
	using BaseComplexT = typename tracking::TrackerTraits<TrackerType>::BaseComplexT;

	using NeededConfigs = detail::TypeList<
								RegenerationConfig,
								TolerancesConfig,
								SharpeningConfig,
								PostProcessingConfig
								>;
};



struct AnyNID : public virtual AnyAlgorithm
{
	virtual ~AnyNID() = default;
};



/**
\brief The Numerical Irreducible Decomposition algorithm.

\note This is framework scaffolding.  The regenerative cascade itself is not yet
implemented; the compute entry points (Run, Solve, RegenerativeCascade) currently
throw.  The class is fully wired for configuration (via detail::Configured and the
NeededConfigs typelist) and observation, and exposes a Tracker and Endgame, so it
slots into the existing Python config interface with no extra plumbing.
*/
template<	typename TrackerType, typename EndgameType,
			typename SystemType >
struct NumericalIrreducibleDecomposition :
					public virtual AnyNID,
					public Observable,
					public detail::Configured<
						typename AlgoTraits< NumericalIrreducibleDecomposition<TrackerType, EndgameType, SystemType> >::NeededConfigs>
{
	// these usings are for getters in python
	using TrackerT = TrackerType;
	using EndgameT = EndgameType;
	using SystemT  = SystemType;


/// a bunch of using statements to reduce typing.
	using BaseComplexT = typename tracking::TrackerTraits<TrackerType>::BaseComplexT;
	using BaseRealT    = typename tracking::TrackerTraits<TrackerType>::BaseRealT;

	using Config = detail::Configured<
						typename AlgoTraits< NumericalIrreducibleDecomposition<TrackerType, EndgameType, SystemType> >::NeededConfigs>;
	using Config::Get;


	using Regeneration   = RegenerationConfig;
	using Tolerances     = TolerancesConfig;
	using Sharpening     = SharpeningConfig;
	using PostProcessing = PostProcessingConfig;

	using ResultT = nag_datatype::NumericalIrreducibleDecomposition<BaseComplexT>;

	// NID owns a cloned target system directly (no policy).
	const SystemType& TargetSystem() const { return target_system_; }
	SystemType&       TargetSystem()       { return target_system_; }


/// constructors

	/**
	Construct a NumericalIrreducibleDecomposition algorithm object from the system to be decomposed.
	*/
	NumericalIrreducibleDecomposition(SystemType const& target)
	 : target_system_(Clone(target)), tracker_(TargetSystem()), endgame_(tracker_)
	{
		DefaultSetup();
	}

	virtual ~NumericalIrreducibleDecomposition() = default;


/// the main functions

	/**
	\brief Main Run() function provided for calling from the blackbox mode.
	*/
	void Run() override
	{
		Solve();
	}

	/**
	\brief Perform the numerical irreducible decomposition.

	\note Not yet implemented -- this is framework scaffolding.
	*/
	void Solve()
	{
		throw std::runtime_error("NumericalIrreducibleDecomposition is not yet implemented");
	}

	/**
	\brief Run the regenerative cascade.

	\note Not yet implemented -- this is framework scaffolding.
	*/
	ResultT RegenerativeCascade()
	{
		throw std::runtime_error("NumericalIrreducibleDecomposition::RegenerativeCascade is not yet implemented");
	}

	/**
	\brief Get the most recently computed decomposition.
	*/
	const ResultT& GetDecomposition() const
	{
		return decomposition_;
	}


/// tracker / endgame access

	const TrackerType& GetTracker() const
	{
		return tracker_;
	}

	TrackerType& GetTracker()
	{
		return tracker_;
	}

	const EndgameType& GetEndgame() const
	{
		return endgame_;
	}

	EndgameType& GetEndgame()
	{
		return endgame_;
	}


/// setup functions

	void DefaultSetup()
	{
		DefaultSettingsSetup();
		DefaultSystemSetup();
	}

	/**
	Fills the configs from default values.
	*/
	void DefaultSettingsSetup()
	{
		this->template Set<Regeneration>(Regeneration());
		this->template Set<Tolerances>(Tolerances());
		this->template Set<Sharpening>(Sharpening());
		this->template Set<PostProcessing>(PostProcessing());
	}

	void DefaultSystemSetup()
	{
		// homogenize + patch the owned target (the old CloneTarget::SystemSetup; a no-op of effect
		// since Solve() throws, but kept so the prepared target is consistent if ever inspected).
		target_system_.Homogenize();
		target_system_.AutoPatch();
	}


private:
	SystemType  target_system_;   ///< the cloned, owned system to be decomposed
	TrackerType tracker_;
	EndgameType endgame_;
	ResultT decomposition_;
};


	} // ns algorithm

} // ns bertini