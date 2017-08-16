//This file is part of Bertini 2.
//
//bertini2/blackbox/global_configs.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/blackbox/global_configs.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/blackbox/global_configs.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/blackbox/global_configs.hpp 

\brief Provides types and utilities for dealing with global defaults for the blackbox routines
*/


#pragma once
#include "bertini2/detail/typelist.hpp"
#include "bertini2/detail/configured.hpp"

#include "bertini2/trackers/config.hpp"
#include "bertini2/endgames/config.hpp"
#include "bertini2/nag_algorithms/common/config.hpp"

namespace bertini{

namespace blackbox{

namespace config {


namespace {
	using namespace tracking;
	using namespace endgame;
	using namespace algorithm;
}

struct Configs
{
	

	using Tracking = detail::TypeList<SteppingConfig, NewtonConfig, FixedPrecisionConfig, AdaptiveMultiplePrecisionConfig, tracking::PrecisionType, Predictor>;

	using Endgame = detail::TypeList<SecurityConfig, EndgameConfig, PowerSeriesConfig, CauchyConfig, TrackBackConfig>;

	template<typename T>
	using Algorithm = detail::TypeList<TolerancesConfig, MidPathConfig, AutoRetrackConfig, SharpeningConfig, RegenerationConfig, PostProcessingConfig, ZeroDimConfig<T>, classic::AlgoChoice>;

	template<typename T>
	using All = detail::ListCat<Tracking, Endgame, Algorithm<T>>;
};


struct Defaults : 
detail::Configured<Configs::All<bertini::dbl>, Configs::All<bertini::mpfr>>
{


};

}
} // namespace blackbox
} // namespace bertini






