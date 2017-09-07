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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire
// Tim Hodges, Colorado State University
// jeb collins, west texas a&m
 
#pragma once

namespace bertini{
	namespace algorithm{

template<typename T>
using SolnCont = std::vector<T>;

namespace classic{

enum class AlgoChoice
{
	EvalFunctions = -4,
	EvalFunctionJacobian = -3,
	NewtonIteration = -2,
	NewtonIterationCondNum = -1,
	ZeroDim = 0,
	NID = 1,
	SampleComponent = 2,
	MembershipTest = 3,
	ExtractWitnessSet = 4,
	WitnessSetProjection = 5,
	IsosingularStab = 6
};

} // namespace classic



struct TolerancesConfig
{	
	using T = NumErrorT;

	T newton_before_endgame = T(1)/T(100000); //E.4.1
	T newton_during_endgame = T(1)/T(1000000); //E.4.2

	T final_tolerance = T(1)/T(100000000000); //E.5.1

	T path_truncation_threshold = T(100000); //E.4.13
};
		
	
struct MidPathConfig
{
	using T = NumErrorT;

	T same_point_tolerance = T(1)/T(100000);
};



struct AutoRetrackConfig
{
	using T = NumErrorT;

	T midpath_decrease_tolerance_factor = T(1)/T(2);
};



struct SharpeningConfig
{
	using T = NumErrorT;

	unsigned sharpendigits; ///< how many digits should be correct after sharpening.
	
	// std::function<Vec<T>> sharpen_method_; ///< function taking a vector, and sharpening it.

	T function_residual_tolerance = Eigen::NumTraits<T>::dummy_precision(); ///< A polynomial (or any function, really) evaluated at a point is considered to be 0 if the magnitude is smaller than this value.  See also RatioTolerance.  **Note that this value depends on the current default precision when this scruct is constructed.**

	T ratio_tolerance = T(99)/T(100); ///<  A computed value is considered to be zero if the ratio of two different approximations is smaller than this value.  See also FunctionTolerance
};



struct RegenerationConfig
{
	using T = NumErrorT;

	bool remove_infinite_endpoints = true; ///<  Bool indicating whether endpoints during the regeneration start point buildup step which are infinite should be discarded.  If you are not interested in infinite solutions, ensure this is true.  RegenRemoveInf

	bool higher_dimension_check = true; ///< RegenHigherDimTest
	unsigned start_level = 0;
	T newton_before_endgame; ///< The tolerance for tracking before reaching the endgame.  SliceTolBeforeEG
	T newton_during_endgame; ///< The tolerance for tracking during the endgame.  SliceTolDuringEG
	T final_tolerance; ///< The final tolerance to track to, using the endgame.  SliceFinalTol
};



struct PostProcessingConfig{
	using T = NumErrorT;
	
	T real_threshold = T(1)/T(100000000); ///< threshold on the imaginary part of a solution being 0.  If the imag part exceeds this, the point is considered complex.  Currently, this is the implemented available way in Bertini2 for determining this, but there are other methods.  Smale's alpha theory provides ways to prove that a point is real.  If this is something you need, please consider adding the method to the library, for all to use!  Or, if this is technically beyond your C++ capabilities, add as an issue on the github page, and indicate it as a feature request.

	T endpoint_finite_threshold = T(1)/T(100000);  ///< The threshold on norm of endpoints being considered infinite.  There is another setting in Tolerances, `path_truncation_threshold`, which tells the path tracker to die if exceeded.  Another related setting is in Security, `max_norm` -- the endgame dies if the norm of the computed approximation exceeds this twice.

	T same_point_tolerance {T(1)/T(10000000000)}; ///< The tolerance for whether two points are the same.  This should be *lower* than the accuracy to which you request your solutions be computed.  Perhaps by at least two orders of magnitude, but the default value is a factor of 10 less stringent.  This also depends on the norm being used to tell whether two points are the same, and the norm used for the convergence condition to terminate tracking.
};

template<typename ComplexT>
struct ZeroDimConfig
{
	unsigned initial_ambient_precision = DoublePrecision();
	unsigned max_num_crossed_path_resolve_attempts = 2; ///< The maximum number of times to attempt to re-solve crossed paths at the endgame boundary.

	ComplexT start_time = ComplexT(1);
	ComplexT endgame_boundary = ComplexT(1)/ComplexT(10);
	ComplexT target_time = ComplexT(0);

	std::string path_variable_name = "ZERO_DIM_PATH_VARIABLE";
};

struct MetaConfig
{
	classic::AlgoChoice tracktype = classic::AlgoChoice::ZeroDim;
};


// a forward declare
template <typename T>
	struct AlgoTraits;


} } // namespaces
