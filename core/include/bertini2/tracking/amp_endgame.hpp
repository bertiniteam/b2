//This file is part of Bertini 2.
//
//amp_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_endgame.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame
// Tim Hodges, Colorado State University



#ifndef BERTINI_TRACKING_AMP_ENDGAME_HPP
#define BERTINI_TRACKING_AMP_ENDGAME_HPP

/**
\file base_endgame.hpp

\brief Contains parent class, Endgame, the parent class for all endgames.
*/


#include <typeinfo>
#include "bertini2/tracking.hpp"
#include "bertini2/system.hpp"
#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include "bertini2/limbo.hpp"
#include <boost/multiprecision/mpfr.hpp>
#include "bertini2/mpfr_complex.hpp"
#include "bertini2/tracking/base_endgame.hpp"

namespace bertini{ 

	namespace tracking {

		namespace endgame {
			

		inline
		void SetPrecision(SampCont<mpfr> & samples, unsigned prec)
		{
			for (auto& s : samples)
				for (unsigned ii=0; ii<s.size(); ii++)
					s(ii).precision(prec);
		}

		inline
		void SetPrecision(TimeCont<mpfr> & times, unsigned prec)
		{
			for (auto& t : times)
				t.precision(prec);
		}

		inline
		unsigned MaxPrecision(SampCont<mpfr> const& samples)
		{
			unsigned max_precision = 0;
			for (auto& s : samples)
				if(Precision(s(0)) > max_precision)
					max_precision = Precision(s(0));
			return max_precision;
		}

		inline
		unsigned MaxPrecision(TimeCont<mpfr> const& times)
		{
			unsigned max_precision = 0;
			for (auto& t : times)
				if(Precision(t) > max_precision)
					max_precision = Precision(t);
			return max_precision;
		}


		//does not a thing, because cannot.
		inline
		unsigned EnsureAtUniformPrecision(TimeCont<dbl> & times, SampCont<dbl> & derivatives)
		{
			return DoublePrecision();
		}


		//changes precision of mpfr to highest needed precision for the samples.
		inline
		unsigned EnsureAtUniformPrecision(TimeCont<mpfr> & times, SampCont<mpfr> & samples)
		{
			unsigned max_precision = max(
			                             mpfr_float::default_precision(),
			                             MaxPrecision(samples)
			                             ); 

			if (max_precision != mpfr_float::default_precision())
				BOOST_LOG_TRIVIAL(severity_level::trace) << "EG changing default precision from " << mpfr_float::default_precision() << " to " << max_precision << std::endl;
			
			mpfr_float::default_precision(max_precision);

			SetPrecision(times, max_precision);
			SetPrecision(samples, max_precision);

			return max_precision;
		}


		//changes precision of mpfr to highest needed precision for the samples.
		inline
		unsigned EnsureAtUniformPrecision(TimeCont<mpfr> & times, SampCont<mpfr> & samples, SampCont<mpfr> & derivatives)
		{
			unsigned max_precision = max(
			                             mpfr_float::default_precision(),
			                             MaxPrecision(samples),
			                             MaxPrecision(derivatives)
			                             ); 

			if (max_precision != mpfr_float::default_precision())
				BOOST_LOG_TRIVIAL(severity_level::trace) << "EG changing default precision from " << mpfr_float::default_precision() << " to " << max_precision << std::endl;
			
			mpfr_float::default_precision(max_precision);

			SetPrecision(times, max_precision);
			SetPrecision(samples, max_precision);
			SetPrecision(derivatives, max_precision);

			return max_precision;
		}

		
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
				using BaseComplexType = typename TrackerTraits<TrackerType>::BaseComplexType;
				using BaseRealType = typename TrackerTraits<TrackerType>::BaseRealType;

				using BCT = BaseComplexType;
				using BRT = BaseRealType;
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




				SuccessCode RefineSample(Vec<mpfr> & result, Vec<mpfr> const& current_sample, mpfr const& current_time)
				{
					using RT = mpfr_float;

					auto refinement_success = this->GetTracker().Refine(result,current_sample,current_time,
					                          	RT(this->Tolerances().final_tolerance)/100,
					                          	this->EndgameSettings().max_num_newton_iterations);

					
					if (refinement_success==SuccessCode::HigherPrecisionNecessary ||
					    refinement_success==SuccessCode::FailedToConverge)
					{
						auto prev_precision = this->Precision();
						auto temp_higher_prec = max(prev_precision,LowestMultiplePrecision())+ PrecisionIncrement();
						mpfr_float::default_precision(temp_higher_prec);
						this->GetTracker().ChangePrecision(temp_higher_prec);


						auto next_sample_higher_prec = current_sample;
						auto result_higher_prec = Vec<mpfr>(current_sample.size());
						auto time_higher_precision = current_time;

						assert(time_higher_precision.precision()==mpfr_float::default_precision());

						refinement_success = this->GetTracker().Refine(result_higher_prec,
						                                               next_sample_higher_prec,
						                                               time_higher_precision,
					                          							RT(this->Tolerances().final_tolerance)/100,
					                          							this->EndgameSettings().max_num_newton_iterations);

						mpfr_float::default_precision(prev_precision);
						this->GetTracker().ChangePrecision(prev_precision);
						result = result_higher_prec;
						assert(result(0).precision()==mpfr_float::default_precision());
					}

					return refinement_success;
				}

				SuccessCode RefineSample(Vec<dbl> & result, Vec<dbl> const& current_sample, dbl const& current_time)
				{
					using RT = double;

					auto refinement_success = this->GetTracker().Refine(result,current_sample,current_time,
					                          	RT(this->Tolerances().final_tolerance)/100,
					                          	this->EndgameSettings().max_num_newton_iterations);

					
					if (refinement_success==SuccessCode::HigherPrecisionNecessary ||
					    refinement_success==SuccessCode::FailedToConverge)
					{
						auto prev_precision = this->Precision();
						auto temp_higher_prec = LowestMultiplePrecision();
						mpfr_float::default_precision(temp_higher_prec);
						this->GetTracker().ChangePrecision(temp_higher_prec);


						auto next_sample_higher_prec = Vec<mpfr>(current_sample.size());
						for (int ii=0; ii<current_sample.size(); ++ii)
							next_sample_higher_prec(ii) = mpfr(current_sample(ii));

						auto result_higher_prec = Vec<mpfr>(current_sample.size());
						mpfr time_higher_precision(current_time);

						refinement_success = this->GetTracker().Refine(result_higher_prec,
						                                               next_sample_higher_prec,
						                                               time_higher_precision,
					                          							RT(this->Tolerances().final_tolerance)/100,
					                          							this->EndgameSettings().max_num_newton_iterations);


						mpfr_float::default_precision(prev_precision);
						this->GetTracker().ChangePrecision(prev_precision);
						for (unsigned ii(0); ii<current_sample.size(); ++ii)
							result(ii) = dbl(result_higher_prec(ii));
					}

					return refinement_success;
				}



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



				explicit AMPEndgame(TrackerType const& tr, const std::tuple< const config::Endgame<BRT>&, const config::Security<BRT>&, const config::Tolerances<BRT>& >& settings )
			      : EndgameBase<TrackerType>(tr, settings),
			        precision_(mpfr_float::default_precision())
			   	{}

			    template< typename... Ts >
   				AMPEndgame(TrackerType const& tr, const Ts&... ts ) : AMPEndgame(tr, Unpermute< config::Endgame<BRT>, config::Security<BRT>, config::Tolerances<BRT> >( ts... ) ) 
   				{}

			};
			
		}// end namespace endgame
	}// end namespace tracking 
}// end namespace bertini
				

#endif
