//This file is part of Bertini 2.0.
//
//explicit_predictors.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//explicit_predictors.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with explicit_predictors.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  copyright 2015-2026
//  James B. Collins, West Texas A&M University, and the Bertini 2 team

/**
 \file explicit_predictors.cpp

 \brief Explicit template instantiations for ExplicitRKPredictor.

 The Butcher tables themselves are now C++17 inline members in the header
 (issue #287); this translation unit exists solely to hold the explicit template
 instantiation definitions that pair with the extern template declarations in the
 header (ADR-0014), keeping each including TU from emitting its own copies.
 */

#include "bertini2/trackers/explicit_predictors.hpp"

namespace bertini{
	namespace tracking{
		namespace predict{

		// Explicit instantiation definitions for the two concrete numeric types.
		// These pair with the extern template declarations in explicit_predictors.hpp
		// and prevent each including TU from emitting its own copy.

		template SuccessCode ExplicitRKPredictor::Predict<complex_dbl>(
		    Vec<complex_dbl>&, StepMetadata&, System const&, Vec<complex_dbl> const&, complex_dbl, complex_dbl const&,
		    unsigned&, unsigned, NumErrorT const&, AdaptiveMultiplePrecisionConfig const*);
		template SuccessCode ExplicitRKPredictor::Predict<complex_mp>(
		    Vec<complex_mp>&, StepMetadata&, System const&, Vec<complex_mp> const&, complex_mp, complex_mp const&,
		    unsigned&, unsigned, NumErrorT const&, AdaptiveMultiplePrecisionConfig const*);

		template SuccessCode ExplicitRKPredictor::FullStep<complex_dbl>(
		    Vec<complex_dbl>&, System const&, Vec<complex_dbl> const&, complex_dbl const&, complex_dbl const&);
		template SuccessCode ExplicitRKPredictor::FullStep<complex_mp>(
		    Vec<complex_mp>&, System const&, Vec<complex_mp> const&, complex_mp const&, complex_mp const&);

		template void ExplicitRKPredictor::SetNormsCond<complex_dbl>(
		    double&, double&, double&, unsigned&, unsigned);
		template void ExplicitRKPredictor::SetNormsCond<complex_mp>(
		    double&, double&, double&, unsigned&, unsigned);

		template SuccessCode ExplicitRKPredictor::SetErrorEstimate<complex_dbl>(double&, complex_dbl const&);
		template SuccessCode ExplicitRKPredictor::SetErrorEstimate<complex_mp>(double&, complex_mp const&);

		template SuccessCode ExplicitRKPredictor::SetSizeProportion<complex_dbl>(double&, complex_dbl const&);
		template SuccessCode ExplicitRKPredictor::SetSizeProportion<complex_mp>(double&, complex_mp const&);

		template SuccessCode ExplicitRKPredictor::EvalRHS<complex_dbl>(
		    System const&, Vec<complex_dbl> const&, complex_dbl const&, Mat<complex_dbl>&, unsigned);
		template SuccessCode ExplicitRKPredictor::EvalRHS<complex_mp>(
		    System const&, Vec<complex_mp> const&, complex_mp const&, Mat<complex_mp>&, unsigned);

		template void ExplicitRKPredictor::FillButcherTable<double>(
		    int, Mat<mpq_rational> const&, Mat<mpq_rational> const&,
		    Mat<mpq_rational> const&, Mat<mpq_rational> const&);
		template void ExplicitRKPredictor::FillButcherTable<real_mp>(
		    int, Mat<mpq_rational> const&, Mat<mpq_rational> const&,
		    Mat<mpq_rational> const&, Mat<mpq_rational> const&);

		template void ExplicitRKPredictor::FillButcherTable<double>(
		    int, Mat<mpq_rational> const&, Mat<mpq_rational> const&, Mat<mpq_rational> const&);
		template void ExplicitRKPredictor::FillButcherTable<real_mp>(
		    int, Mat<mpq_rational> const&, Mat<mpq_rational> const&, Mat<mpq_rational> const&);

		} // re: predict
	}// re: tracking
}// re: bertini





