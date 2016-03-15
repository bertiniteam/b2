//This file is part of Bertini 2.0.
//
//base_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//base_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with predict.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  base_endgame.hpp
//
//  copyright 2015
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//  Fall 2015


#ifndef BERTINI_TRACKING_AMP_ENDGAME_HPP
#define BERTINI_TRACKING_AMP_ENDGAME_HPP

/**
\file base_endgame.hpp

\brief Contains parent class, Endgame, the parent class for all endgames.
*/


#include <typeinfo>
#include "tracking.hpp"
#include "system.hpp"
#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include "limbo.hpp"
#include <boost/multiprecision/mpfr.hpp>
#include "mpfr_complex.hpp"
#include "tracking/base_endgame.hpp"

namespace bertini{ 

	namespace tracking {

		namespace endgame {
			


		

		/**
		\class Endgame

		\brief Base endgame class for all endgames offered in Bertini2.
	
		\see PowerSeriesEndgame
		\see CauchyEndgame
		
		## Using an endgame

		Endgames in Bertini2 are the engine for finishing homotopy continuation where we may encounter singular solutions.
		The path is implicitly described by the system being tracked.

		## Purpose 

		Since the Bertini Endgames have common functionality, and we want to be able to call arbitrary algorithms using and tracker type, we use inheritance. That is, there is common functionality in all endgames, such as

		ComputeInitialSamples

		Also, there are settings that will be kept at this level to not duplicate code. 
			
		## Creating a new endgame type

		 To create a new endgame type, inherit from this class. 
		*/
			
			template<class TrackerType>
			class AMPEndgame : public EndgameBase<TrackerType>
			{

			protected:
				

				mutable unsigned precision_;
				unsigned initial_precision_;
				bool preserve_precision_ = true; ///< Whether the endgame should change back to the initial precision after running 


				/**
				\brief Change precision of all temporary internal state variables.

				This excludes those which cannot be re-written without copying.

				\brief new_precision The new precision to adjust to.
				*/
				void ChangeTemporariesPrecision(unsigned new_precision) const
				{
					ChangeTemporariesPrecisionImpl(new_precision);
				}

				/**
				\brief Converts from multiple to different precision multiple precision

				Changes the precision of the internal temporaries to desired precision

				\param new_precision The new precision.
				*/
				void MultipleToMultiple(unsigned new_precision) const
				{
					mpfr_float::default_precision(new_precision);
					precision_ = new_precision;

					const auto& source_point = std::get<Vec<mpfr> >(this->final_approximation_at_origin_);
					auto& target_point = std::get<Vec<mpfr> >(this->final_approximation_at_origin_);
					target_point.resize(source_point.size());

					for (unsigned ii=0; ii<source_point.size(); ii++)
					{
						target_point(ii).precision(new_precision);
						target_point(ii) = mpfr(source_point(ii));
					}
					ChangeTemporariesPrecision(new_precision);

					MultipleToMultipleImpl(new_precision);
				}

				/**
				\brief Converts from multiple to different precision multiple precision

				Changes the precision of the internal temporaries to desired precision

				\param new_precision The new precision.
				*/
				void DoubleToMultiple(unsigned new_precision) const
				{
					mpfr_float::default_precision(new_precision);
					precision_ = new_precision;
					const auto& source_point = std::get<Vec<dbl> >(this->final_approximation_at_origin_);
					auto& target_point = std::get<Vec<mpfr> >(this->final_approximation_at_origin_);

					target_point.resize(source_point.size());
					for (unsigned ii=0; ii<source_point.size(); ii++)
					{
						target_point(ii).precision(new_precision);
						target_point(ii) = mpfr(source_point(ii));
					}
					ChangeTemporariesPrecision(new_precision);

					DoubleToMultipleImpl(new_precision);
				}

				/**
				\brief Converts from multiple to different precision multiple precision

				Changes the precision of the internal temporaries to desired precision

				\param new_precision The new precision.
				*/
				void MultipleToDouble() const
				{
					mpfr_float::default_precision(DoublePrecision());
					precision_ = DoublePrecision();

					const auto& source_point = std::get<Vec<mpfr> >(this->final_approximation_at_origin_);
					auto& target_point = std::get<Vec<dbl> >(this->final_approximation_at_origin_);
					target_point.resize(source_point.size());

					for (unsigned ii=0; ii<source_point.size(); ii++)
						target_point(ii) = dbl(source_point(ii));

					ChangeTemporariesPrecision(DoublePrecision());

					MultipleToDoubleImpl();
				}

				bool PrecisionSanityCheck() const
				{
					return true;
				}


				virtual void ChangeTemporariesPrecisionImpl(unsigned new_precision) const = 0;
				virtual void MultipleToMultipleImpl(unsigned new_precision) const = 0;
				virtual void DoubleToMultipleImpl(unsigned new_precision) const = 0;
				virtual void MultipleToDoubleImpl() const = 0;

			public:	

				/**
				Change precision of endgame object to next_precision.  Converts the internal temporaries, and adjusts precision of system. 

				\param new_precision The precision to change to.
				\return SuccessCode indicating whether the change was successful.  If the precision increases, and the refinement loop fails, this could be not Success.  Changing down is guaranteed to succeed.
				*/
				SuccessCode ChangePrecision(unsigned new_precision) const
				{
					if (new_precision==precision_) // no op
						return SuccessCode::Success;

					if (new_precision==DoublePrecision() && precision_>DoublePrecision())
						MultipleToDouble();
					else if(new_precision > DoublePrecision() && precision_ == DoublePrecision())
						DoubleToMultiple(new_precision);
					else
						MultipleToMultiple(new_precision);


					#ifndef BERTINI_DISABLE_ASSERTS
					assert(PrecisionSanityCheck() && "precision sanity check failed.  some internal variable is not in correct precision");
					#endif

					return SuccessCode::Success;
				}

				auto Precision() const
				{ return precision_; }



				explicit AMPEndgame(TrackerType const& tr, const std::tuple< const config::Endgame&, const config::Security&, const config::Tolerances& >& settings )
			      : EndgameBase<TrackerType>(tr, settings),
			        precision_(mpfr_float::default_precision())
			   	{}

			    template< typename... Ts >
   				AMPEndgame(TrackerType const& tr, const Ts&... ts ) : AMPEndgame(tr, Unpermute< config::Endgame, config::Security, config::Tolerances >( ts... ) ) 
   				{}

			};
			
		}// end namespace endgame
	}// end namespace tracking 
}// end namespace bertini
				

#endif
