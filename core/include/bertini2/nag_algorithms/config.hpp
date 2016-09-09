//This file is part of Bertini 2.
//
//nag_algorithms/config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//nag_algorithms/config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with nag_algorithms/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame
// Tim Hodges, Colorado State University
 
#pragma once

namespace bertini{
	namespace algorithm{

template<typename T>
using SolnCont = std::vector<T>;

		namespace config{


template<typename T>
struct Tolerances
{	
	T newton_before_endgame = T(1)/T(100000); //E.4.1
	T newton_during_endgame = T(1)/T(1000000); //E.4.2

	T final_tolerance = T(1)/T(100000000000); //E.5.1
	T final_tolerance_multiplier = T(10); // This multiplier is used to cluster or de-cluster points at the target system. 

	T path_truncation_threshold = T(100000); //E.4.13
	// T final_tolerance_times_final_tolerance_multiplier = final_tolerance * final_tolerance_multiplier;
};
			

			template<typename T>
			struct AutoRetrack
			{
				T midpath_decrease_tolerance_factor = T(1)/T(2);
				T boundary_near_tol = T(1)/T(100000);
			};


template<typename T>
struct Sharpening
{
	unsigned sharpendigits; ///< how many digits should be correct after sharpening.
	
	std::function<Vec<T>> sharpen_method_; ///< function taking a vector, and sharpening it.

	T function_residual_tolerance = Eigen::NumTraits<T>::dummy_precision(); ///< A polynomial (or any function, really) evaluated at a point is considered to be 0 if the magnitude is smaller than this value.  See also RatioTolerance

	T ratio_tolerance = T(99)/T(100); ///<  A computed value is considered to be zero if the ratio of two different approximations is smaller than this value.  See also FunctionTolerance
};


template<typename T>
struct Regeneration
{
	bool remove_infinite_endpoints = true; ///<  Bool indicating whether endpoints during the regeneration start point buildup step which are infinite should be discarded.  If you are not interested in infinite solutions, ensure this is true.  RegenRemoveInf

	bool higher_dimension_check = true; ///< RegenHigherDimTest

	T newton_before_endgame; ///< The tolerance for tracking before reaching the endgame.  SliceTolBeforeEG
	T newton_during_endgame; ///< The tolerance for tracking during the endgame.  SliceTolDuringEG
	T final_tolerance; ///< The final tolerance to track to, using the endgame.  SliceFinalTol
};


template<typename T>
struct PostProcessing{
	T real_threshold = T(1)/T(100000000); ///< threshold on the imaginary part of a solution being 0.  If the imag part exceeds this, the point is considered complex.  Currently, this is the implemented available way in Bertini2 for determining this, but there are other methods.  Smale's alpha theory provides ways to prove that a point is real.  If this is something you need, please consider adding the method to the library, for all to use!  Or, if this is technically beyond your C++ capabilities, add as an issue on the github page, and indicate it as a feature request.

	T endpoint_finite_threshold = T(1)/T(100000);  ///< The threshold on norm of endpoints being considered infinite.  There is another setting in Tolerances, `path_truncation_threshold`, which tells the path tracker to die if exceeded.  Another related setting is in Security, `max_norm` -- the endgame dies if the norm of the computed approximation exceeds this twice.

	T final_tol_multiplier{10}; ///< Multiply this value by FinalTol to yield the tolerance used to numerically determine if two endpoints should be considered the same point.
};


struct ZeroDim
{
	unsigned max_num_crossed_path_resolve_attempts = 2; ///< The maximum number of times to attempt to re-solve crossed paths at the endgame boundary.
};


		}
	}
}
