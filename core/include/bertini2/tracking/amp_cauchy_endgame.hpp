//This file is part of Bertini 2.0.
//
//cauchy_endgame.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//cauchy_endgame.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with predict.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  cauchy_endgame.hpp
//
//  copyright 2015
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//  Spring 2015


#ifndef BERTINI_TRACKING_CAUCHY_HPP
#define BERTINI_TRACKING_CAUCHY_HPP

/**
\file cauchy_endgame.hpp

\brief Contains the cauchy endgame type, and necessary functions.
*/


#include <boost/multiprecision/gmp.hpp>
#include <iostream>
#include "tracking/amp_endgame.hpp"
#include <cstdio>
#include <typeinfo> 
#include <deque>
#include <algorithm> 



namespace bertini
{ 

	namespace tracking 
	{

		namespace endgame 
		{

		/** 
		\class CauchyEndgame

		\brief Class used to finish tracking paths during Homotopy Continuation.
		
		
		## Explanation
		
		The bertini::CauchyEngame class enables us to finish tracking on possibly singular paths on an arbitrary square homotopy.  

		The intended usage is to:

		1. Create a system, tracker, and instantiate some settings.
		2. Using the tracker created track to the engame boundary. 
		3. Create a CauchyEndgame, associating it to the system you are going to solve or track on.
		4. For each path being tracked send the CauchyEndgame the time value and other variable values at that time. 
		5. The CauchyEndgame, if successful, will store the target systems solution at $t = 0$.

		## Example Usage
		Below we demonstrate a basic usage of the CauchyEndgame class to find the singularity at $t = 0$. 

		The pattern is as described above: create an instance of the class, feeding it the system to be used, and the endgame boundary time and other variable values at the endgame boundary. 

		\code{.cpp}
		mpfr_float::default_precision(30); // set initial precision.  This is not strictly necessary.
		using namespace bertini::tracking;

		// 1. Create the system
		System sys;
		Var x = std::make_shared<Variable>("x"), t = std::make_shared<Variable>("t"), y = std::make_shared<Variable>("y");
		VariableGroup vars{x,y};
		sys.AddVariableGroup(vars); 
		sys.AddPathVariable(t);
		// Define homotopy system
		sys.AddFunction((pow(x-1,3))*(1-t) + (pow(x,3) + 1)*t);
		sys.AddFunction((pow(y-1,2))*(1-t) + (pow(y,2) + 1)*t);

		//2. Setup a tracker. 
		bertini::tracking::AMPTracker tracker(sys);
	
		bertini::tracking::config::Stepping stepping_preferences;
		bertini::tracking::config::Newton newton_preferences;

		tracker.Setup(bertini::tracking::config::Predictor::Euler,
                mpfr_float("1e-5"),
                mpfr_float("1e5"),
                stepping_preferences,
                newton_preferences);
	
		tracker.AMPSetup(AMP);

		// We make the assumption that we have tracked to t = 0.1 (value at t = 0.1 known from Bertini 1.5)
		mpfr current_time(1);
		Vec<mpfr> current_space(2);
		current_time = mpfr(".1");
		current_space <<  mpfr("5.000000000000001e-01", "9.084258952712920e-17") ,mpfr("9.000000000000001e-01","4.358898943540673e-01");

		//  3. (Optional) configure settings for specific endgame. 
		bertini::tracking::config::Endgame endgame_settings;
		bertini::tracking::config::Cauchy cauchy_settings;
		bertini::tracking::config::Security endgame_security_settings;

		//4. Create the CauchyEngame object, by sending in the tracker and any other settings. Notice the endgame is templated by the tracker. 
		bertini::tracking::endgame::CauchyEndgame<bertini::tracking::AMPTracker> My_Endgame(tracker,cauchy_settings,endgame_settings,endgame_security_settings);

		//Calling the CauchyEG member function actually runs the endgame. 
		My_Endgame.CauchyEG(current_time,current_space);

		//Access solution at t = 0, using the .get_final_approximation_at_origin() member function. 
		auto Answer  = My_Endgame.get_final_approximation_at_origin();
		\endcode
		
		If this documentation is insufficient, please contact the authors with suggestions, or get involved!  Pull requests welcomed.
		
		## Testing
		Test suite driving this class: endgames_test.

		File: test/endgames/cauchy_class_test.cpp

		Functionality tested: All member functions of the CauchyEndgame have been tested in variable precision. There is also, test running different systems to find singular solutions at t = 0.


		*/



		
			template<typename TrackerType> 
			class CauchyEndgame : public AMPEndgame<TrackerType>
			{
			private:

				/**
				\brief Settings that are specific to the Cauchy endgame. 
				*/	
				config::Cauchy cauchy_settings_; 

				/**
				\brief A deque of times that are specifically used to compute the power series approximation for the Cauchy endgame. 
				*/
				mutable std::tuple<TimeCont<dbl>, TimeCont<mpfr> > pseg_times_;
				/**
				\brief A deque of samples that are in correspondence with the pseg_times_. These samples are also used to compute the first power series approximation for the Cauchy endgame.
				*/
				mutable std::tuple<SampCont<dbl>, SampCont<mpfr> > pseg_samples_; //samples used for the first approximation. 
				/**
				\brief A deque of times that are collected by CircleTrack. These samples are used to compute all approximations of the origin after the first. 
				*/
				mutable std::tuple<TimeCont<dbl>, TimeCont<mpfr> > cauchy_times_;
				/**
				\brief A deque of samples collected by CircleTrack. Computed a mean of the values of this deque, after a loop has been closed, will give an approximation of the origin.
				*/
				mutable std::tuple<SampCont<dbl>, SampCont<mpfr> > cauchy_samples_;



				
				void ChangeTemporariesPrecisionImpl(unsigned new_precision) const override
				{ }

				void MultipleToMultipleImpl(unsigned new_precision) const override
				{
					{
						auto& times = std::get<TimeCont<mpfr> >(cauchy_times_);
						for (auto& t : times)
							t.precision(new_precision);

						auto& samples = std::get<SampCont<mpfr> >(cauchy_samples_);
						for (auto& s : samples)
							for (unsigned ii=0; ii<s.size(); ++ii)
								s(ii).precision(new_precision);
					}
					{
						auto& times = std::get<TimeCont<mpfr> >(pseg_times_);
						for (auto& t : times)
							t.precision(new_precision);

						auto& samples = std::get<SampCont<mpfr> >(pseg_samples_);
						for (auto& s : samples)
							for (unsigned ii=0; ii<s.size(); ++ii)
								s(ii).precision(new_precision);
					}
				}

				void DoubleToMultipleImpl(unsigned new_precision) const override
				{
					{
						auto& times_m = std::get<TimeCont<mpfr> >(cauchy_times_);
						auto& times_d = std::get<TimeCont<dbl > >(cauchy_times_);

						for (unsigned ii=0; ii<times_m.size(); ++ii)
						{
							times_m[ii].precision(new_precision);
							times_m[ii] = mpfr(times_d[ii]);
						}


						auto& samples_m = std::get<SampCont<mpfr> >(cauchy_samples_);
						auto& samples_d = std::get<SampCont<dbl > >(cauchy_samples_);

						for (unsigned ii=0; ii<times_m.size(); ++ii)
						{
							auto& s = samples_m[ii];
							auto& sd = samples_d[ii];
							for (unsigned jj=0; jj<s.size(); ++jj)
							{
								s(jj).precision(new_precision);
								s(jj) = mpfr(sd(jj));
							}
						}
					}
					{
						auto& times_m = std::get<TimeCont<mpfr> >(pseg_times_);
						auto& times_d = std::get<TimeCont<dbl> >(pseg_times_);

						for (unsigned ii=0; ii<times_m.size(); ++ii)
						{
							times_m[ii].precision(new_precision);
							times_m[ii] = mpfr(times_d[ii]);
						}


						auto& samples_m = std::get<SampCont<mpfr> >(pseg_samples_);
						auto& samples_d = std::get<SampCont<dbl > >(pseg_samples_);

						for (unsigned ii=0; ii<times_m.size(); ++ii)
						{
							auto& s = samples_m[ii];
							auto& sd = samples_d[ii];
							for (unsigned jj=0; jj<s.size(); ++jj)
							{
								s(jj).precision(new_precision);
								s(jj) = mpfr(sd(jj));
							}
						}
					}
				}

				void MultipleToDoubleImpl() const override
				{
					{
						auto& times_m = std::get<TimeCont<mpfr> >(cauchy_times_);
						auto& times_d = std::get<TimeCont<dbl> >(cauchy_times_);

						for (unsigned ii=0; ii<times_m.size(); ++ii)
						{
							times_d[ii] = dbl(times_m[ii]);
						}


						auto& samples_m = std::get<SampCont<mpfr> >(cauchy_samples_);
						auto& samples_d = std::get<SampCont<dbl > >(cauchy_samples_);

						for (unsigned ii=0; ii<times_m.size(); ++ii)
						{
							auto& s = samples_m[ii];
							auto& sd = samples_d[ii];
							for (unsigned jj=0; jj<s.size(); ++jj)
							{
								sd(jj) = dbl(s(jj));
							}
						}
					}
					{
						auto& times_m = std::get<TimeCont<mpfr> >(pseg_times_);
						auto& times_d = std::get<TimeCont<dbl> >(pseg_times_);

						for (unsigned ii=0; ii<times_m.size(); ++ii)
						{
							times_d[ii] = dbl(times_m[ii]);
						}


						auto& samples_m = std::get<SampCont<mpfr> >(pseg_samples_);
						auto& samples_d = std::get<SampCont<dbl > >(pseg_samples_);

						for (unsigned ii=0; ii<times_m.size(); ++ii)
						{
							auto& s = samples_m[ii];
							auto& sd = samples_d[ii];
							for (unsigned jj=0; jj<s.size(); ++jj)
							{
								sd(jj) = dbl(s(jj));
							}
						}
					}
				}

			public:


				// // some getters
				// const std::deque< Vec<mpfr> >& CauchySamples() const {return cauchy_samples_;}
				// const std::deque<mpfr>& CauchyTimes() const {return cauchy_times_;}

				/**
				\brief Function that clears all samples and times from data members for the Cauchy endgame
				*/
				template<typename CT>
				void ClearTimesAndSamples()
				{
					std::get<TimeCont<CT> >(pseg_times_).clear(); 
					std::get<TimeCont<CT> >(cauchy_times_).clear(); 
					std::get<SampCont<CT> >(pseg_samples_).clear(); 
					std::get<SampCont<CT> >(cauchy_samples_).clear();}
				/**
				\brief Setter for the deque holding time values for the power series approximation of the Cauchy endgame. 
				*/
				template<typename CT>
				void SetPSEGTimes(TimeCont<CT> pseg_times_to_set) 
				{ std::get<TimeCont<CT> >(pseg_times_) = pseg_times_to_set;}

				/**
				\brief Getter for the deque holding time values for the power series approximation of the Cauchy endgame.
				*/
				template<typename CT>
				TimeCont<CT> GetPSEGTimes() {return std::get<TimeCont<CT> >(pseg_times_);}

				/**
				\brief Setter for the deque holding space values for the power series approximation of the Cauchy endgame. 
				*/
				template<typename CT>
				void SetPSEGSamples(SampCont<CT> pseg_samples_to_set) { std::get<SampCont<CT> >(pseg_samples_) = pseg_samples_to_set;}

				/**
				\brief Getter for the deque holding space values for the power series approximation of the Cauchy endgame. 
				*/
				template<typename CT>
				SampCont<CT> GetPSEGSamples() {return std::get<SampCont<CT> >(pseg_samples_);}

				/**
				\brief Setter for the deque holding space values for the Cauchy endgame. 
				*/
				template<typename CT>
				void SetCauchySamples(SampCont<CT> cauchy_samples_to_set) { std::get<SampCont<CT> >(cauchy_samples_) = cauchy_samples_to_set;}

				/**
				\brief Getter for the deque holding space values for the Cauchy endgame. 
				*/
				template<typename CT>
				SampCont<CT> GetCauchySamples() {return std::get<SampCont<CT> >(cauchy_samples_);}

				/**
				\brief Setter for the deque holding time values for the Cauchy endgame. 
				*/
				template<typename CT>
				void SetCauchyTimes(TimeCont<CT> cauchy_times_to_set) { std::get<TimeCont<CT> >(cauchy_times_) = cauchy_times_to_set;}

				/**
				\brief Getter for the deque holding time values for the Cauchy endgame. 
				*/
				template<typename CT>
				TimeCont<CT> GetCauchyTimes() {return std::get<TimeCont<CT> >(cauchy_times_);}
	
				
				/**
				\brief Setter for the specific settings in tracking_conifg.hpp under Cauchy.
				*/
				void SetCauchySettings(config::Cauchy new_cauchy_settings){cauchy_settings_ = new_cauchy_settings;}

				/**
				\brief Getter for the specific settings in tracking_conifg.hpp under Cauchy.
				*/
				config::Cauchy GetCauchySettings(){return cauchy_settings_;}
				

				explicit CauchyEndgame(TrackerType const& tr, 
				                            const std::tuple< config::Cauchy const&, 
				                            					const config::Endgame&, 
				                            					const config::Security&, 
				                            					const config::Tolerances& >& settings )
			      : AMPEndgame<TrackerType>(tr, std::get<1>(settings), std::get<2>(settings), std::get<3>(settings) ), 
			          cauchy_settings_( std::get<0>(settings) )
			   	{ }

			    template< typename... Ts >
   				CauchyEndgame(TrackerType const& tr, const Ts&... ts ) : CauchyEndgame(tr, Unpermute<config::Cauchy, config::Endgame, config::Security, config::Tolerances >( ts... ) ) 
   				{}


				~CauchyEndgame() {};


				/*
				Input: 
					A starting time and starting sample space to start tracking. 

				Output: A successcode to see if we have successfully tracked around the origin or if we have encountered an error. 

				Details:
					This function uses the number of sample points to track in a polygonal path around the origin. 
					This function should be called the number of times it takes to loop around the origin and get back the original 
					point we started with. In essence this function will help determine the cycle number. A deque of cauchy_samples and cauchy_times is populated in this function. 
				*/
				/**
				\brief A function that will track around the origin once. 

				\param[out] starting_time The time value at which we start trackin around the origin. 
				\param starting_sample The current variable values at staring_time.

				\tparam CT The complex number type.
				*/

				template<typename CT> 
				SuccessCode CircleTrack(CT const& starting_time, Vec<CT> const& starting_sample)
				{	
					using RT = typename Eigen::NumTraits<CT>::Real;
					using std::acos;

					if (this->EndgameSettings().num_sample_points < 3) // need to make sure we won't track right through the origin.
					{
						std::stringstream err_msg;
						err_msg << "ERROR: The number of sample points " << this->EndgameSettings().num_sample_points << " for circle tracking must be >= 3";
						throw std::runtime_error(err_msg.str());
					}	

					// the initial sample has already been added to the sample repo... so don't do that here, please
					
					RT radius = abs(starting_time), angle = arg(starting_time);
					RT angle_increment = 2*acos(RT(-1)) / (this->EndgameSettings().num_sample_points);

					auto number_of_steps = 0;
					CT i(0,1); 

					Vec<CT> current_sample = starting_sample;
					Vec<CT> next_sample;
					CT current_time = starting_time;
					CT next_time;

				

					for (unsigned ii = 0; ii < this->EndgameSettings().num_sample_points; ++ii)
					{
						//setting up the time value for the next sample. 
						if (ii!=this->EndgameSettings().num_sample_points-1)
						{
							using std::polar;
							using bertini::polar;
							angle += angle_increment;
							next_time = polar(radius, angle);
						}
						else
							next_time = starting_time;

						auto tracking_success = this->GetTracker().TrackPath(next_sample, current_time, next_time, current_sample);	
						if (tracking_success != SuccessCode::Success)
						{
							std::cout << "tracker fail in circle track, type " << int(tracking_success) << std::endl;
							return tracking_success;
						}
						
						auto refinement_success = this->RefineSample(next_sample, next_sample, next_time);
						if (refinement_success != SuccessCode::Success)
						{
							std::cout << "refinement fail in circle track, type " << int(refinement_success) << std::endl;
							return refinement_success;
						}
						

						current_sample = next_sample;
						current_time = next_time;
						std::get<TimeCont<CT> >(cauchy_times_).push_back(current_time);
						std::get<SampCont<CT> >(cauchy_samples_).push_back(current_sample);
					}

					return SuccessCode::Success;

				}//end CircleTrack


				
				/*
				Input: 
					All input is available through class data members. 

				Output: An mpfr_float that is an estimate on the c/k ratio shown in the book under Cauchy Endgame.

				Details:
					Using the formula for the cycle test outline in the Bertini Book pg. 53, we can compute an estimate for c/k.
					This estimate is used to help find stabilization of the cycle number. 
				*/

				/**
				\brief A function that uses the assumption of being in the endgame operating zone to compute an approximation of the ratio c over k. 
				When the cycle number stabilizes we will see that the different approximations of c over k will stabilize. 
				Returns the computed value of c over k. 
				*/
				template<typename CT>
				auto ComputeCOverK() -> typename Eigen::NumTraits<CT>::Real
				{//Obtain samples for computing C over K.
					using RT = typename Eigen::NumTraits<CT>::Real;
					using std::abs;
					using std::log;

					const auto& pseg_samples = std::get<SampCont<CT> >(pseg_samples_);

					assert(pseg_samples.size()>=3);
					const Vec<CT> & sample0 = pseg_samples[0];
					const Vec<CT> & sample1 = pseg_samples[1];
					const Vec<CT> & sample2 = pseg_samples[2];

					Vec<CT> rand_vector = Vec<CT>::Random(sample0.size()); //should be a row vector for ease in multiplying.


					// //DO NOT USE Eigen .dot() it will do conjugate transpose which is not what we want.
					// //Also, the .transpose*rand_vector returns an expression template that we do .norm of since abs is not available for that expression type. 
					RT estimate = abs(log(abs((((sample2 - sample1).transpose()*rand_vector).norm())/(((sample1 - sample0).transpose()*rand_vector).norm()))));
					estimate = abs(log(RT(this->EndgameSettings().sample_factor)))/estimate;
					if (estimate < 1)
					  	return RT(1);
					else
						return estimate;

				}//end ComputeCOverK

				/*
				Input: An array of c/k estimates. If we have stabilization than these estimates will be withing some treshold. 

				Output: Boolean of true or false. 

				Details: We consider the ratio of consecutive c/k estimates. If we are less than the minimum needed for stabilization we return false.
				Otherwise, we return true. 
				*/

				/**
				\brief Function to determine if ratios of c/k estimates are withing a user defined threshold. 			

				\param c_over_k_array A deque of computed values for c/k.
				*/
				template<typename RT>
				bool CheckForCOverKStabilization(TimeCont<RT> const& c_over_k_array)
				{	
					using std::abs;

					assert(c_over_k_array.size()>=cauchy_settings_.num_needed_for_stabilization);
					for(unsigned ii = 1; ii < cauchy_settings_.num_needed_for_stabilization ; ++ii)
					{
						auto a = abs(c_over_k_array[ii-1]);
						auto b = abs(c_over_k_array[ii]);

						auto divide = a;

						if(a < b)
							divide = a/b;
						else
							divide = b/a;

						if(divide <  cauchy_settings_.minimum_for_c_over_k_stabilization)
							return false;
					}
					return true;

				}//end CheckForCOverKStabilization


				/*
				Input: A time value and the space value above that time.
					
				Output: An mpfr_float representing a tolerance threshold for declaring a loop to be closed. 


				Details: Used in Bertini 1 as a heuristic for computing separatedness of roots. Decided to not be used since assumptions for this tolerance are not usually met. 


				template<typename CT>
				mpfr_float FindToleranceForClosedLoop(CT x_time, Vec<CT> x_sample)
				{
					auto degree_max = std::max(this->GetTracker().AMP_config_.degree_bound,mpfr_float("2.0")); 
					auto K = this->GetTracker().AMP_config_.coefficient_bound;
					mpfr_float N;
					mpfr_float M;
					mpfr_float L;

					if(max_closed_loop_tolerance_ < min_closed_loop_tolerance_)
					{
						max_closed_loop_tolerance_ = min_closed_loop_tolerance_;
					}

					auto error_tolerance = mpfr_float("1e-13");
					if(x_sample.size() <= 1)
					{
						N = degree_max;
					}
					else
					{
						N = ComputeCombination(degree_max + x_sample[0].precision() - 1, x_sample[0].precision() - 1);
					}
					M = degree_max * (degree_max - 1) * N;
					auto jacobian_at_current_time = this->GetSystem().Jacobian(x_sample,x_time);


					auto minimum_singular_value = Eigen::JacobiSVD< Mat<CT> >(jacobian_at_current_time).singularValues()(this->GetSystem().NumVariables() - 1 );

					auto norm_of_sample = x_sample.norm();
					L = pow(norm_of_sample,degree_max - 2);
					auto tol = K * L * M;
					// std::cout << "TOL IS " << tol << '\n';
					if (tol == 0) // fail-safe
   							tol = minimum_singular_value;
 					else
   					{
   					tol = mpfr_float("2.0") / tol;
    				tol = tol * minimum_singular_value;
  					}
 					// make sure that tol is between min_tol & max_tol
					if (tol > max_closed_loop_tolerance_) // tol needs to be <= max_tol
 				    tol = max_closed_loop_tolerance_;
					if (tol < min_closed_loop_tolerance_) // tol needs to be >= min_tol
  					tol = min_closed_loop_tolerance_;
 					return tol;
				}// end FindToleranceForClosedLoop
				*/


				/*
				Input: No input, tolerance used is tracking tolerance during endgame. 
					
				Output: A boolean declaring if we have closed the loop or not. 

				Details: Take a point used to start CircleTrack along with last point tracked to in CircleTrack. We see if the difference between the two points
				is less than the tracking tolerance during endgame. If so, we declare these points to be the same. 	
				*/

				/**
				\brief Function that determines if we have closed a loop after calling CircleTrack().
				*/
				template<typename CT>
				bool CheckClosedLoop()
				{	
					using RT = typename Eigen::NumTraits<CT>::Real;
					auto& times = std::get<TimeCont<CT> >(cauchy_times_);
					auto& samples = std::get<SampCont<CT> >(cauchy_samples_);

					if((samples.front() - samples.back()).norm() < this->GetTracker().TrackingTolerance())
					{
						return true;
					}

					this->GetTracker().Refine(samples.front(),samples.front(),times.front(),RT(this->Tolerances().final_tolerance),this->EndgameSettings().max_num_newton_iterations);
					this->GetTracker().Refine(samples.back(),samples.back(),times.back(),RT(this->Tolerances().final_tolerance),this->EndgameSettings().max_num_newton_iterations);

					if((samples.front() - samples.back()).norm() < this->GetTracker().TrackingTolerance())
					{
						return true;
					}
					return false;	

				}//end CheckClosedLoop

				/*
				Input: No input, all data needed are class data members. 

				Output: Boolean declaring if our ratios are close enough or not. 

				Details: Finds minimum and maximum norms from tracking around the origin. Then we check against the user defined setting minimum_c_over_k_for_stabilization to see if 
				we have good ratios or not. 
					 
				*/

				/**
				\brief After we have used CircleTrack and have successfully closed the loop using CheckClosedLoop we need to check the maximum and minimum norms of the samples collected. 

				If the ratio of the maximum and minimum norm are within the threshold maximum_cauchy_ratio, and the difference is greater than final tolerance than we are successful. 
				*/
				template<typename CT>
				bool RatioEGOperatingZoneTest()
				{	
					using RT = typename Eigen::NumTraits<CT>::Real;
					RT min(1e300);
					RT max(0);
					auto& times = std::get<TimeCont<CT> >(cauchy_times_);
					auto& samples = std::get<SampCont<CT> >(cauchy_samples_);
					if(norm(times.front()) < cauchy_settings_.ratio_cutoff_time)
					{
						return true;
					}
					else
					{
						RT norm;
						for(unsigned int ii=0; ii < this->EndgameSettings().num_sample_points; ++ii)
						{
							norm = samples[ii].norm();
							if(norm > max)
							{
								max = norm;
							}
							if(norm < min)
							{
								min = norm;
							}
						}

						if(min > this->Tolerances().final_tolerance && max > this->Tolerances().final_tolerance)
						{
							norm = min / max;
							if(norm < cauchy_settings_.maximum_cauchy_ratio && (max - min) > this->Tolerances().final_tolerance)
							{
								return false;
							}
						}
					}
					return true;

				}//end RatioEGOperatingZoneTest


				/*
				Input: All input needed is available as class data members. 

				Output: SuccessCode declaring if we are successful or not. 

				Details: Computes cauchy samples using Circle Track. Compares the ratios of the maximum and minimum norms of cauchy samples 
						 using RatioEGOperatingZoneTest. Then attempts to see if we have closed our loop, otherwise it continues the process outlined.
				*/

				/**
				\brief This function tracks into origin while computing loops around the origin. The function is checking to make sure we have reached a time the ratio of the 
				maximum and minimum norms are withing some tolerance. When this is the case we return success. If this does not happen we will return an error depending on the error 
				encountered. 
				*/		 
				template<typename CT>
				SuccessCode InitialCauchyLoops()
				{	
					using RT = typename Eigen::NumTraits<CT>::Real;
					auto& cau_times = std::get<TimeCont<CT> >(cauchy_times_);
					auto& cau_samples = std::get<SampCont<CT> >(cauchy_samples_);
					auto& ps_times = std::get<TimeCont<CT> >(pseg_times_);
					auto& ps_samples = std::get<SampCont<CT> >(pseg_samples_);

					using std::max;
					bool continue_loop = true;

					auto fail_safe_max_cycle_number = max(cauchy_settings_.fail_safe_maximum_cycle_number,this->CycleNumber());

					auto return_value = SuccessCode::Success;

					while (continue_loop)
					{	
						this->CycleNumber(0);
						cau_times.clear();
						cau_samples.clear();
						cau_samples.push_back(ps_samples.back()); // cauchy samples and times should be empty before this point. 
						cau_times.push_back(ps_times.back());

						auto next_time = ps_times.back();
						auto next_sample = ps_samples.back();

						// track around a circle once.  we'll use it to measure whether we believe we are in the eg operating zone, based on the ratio of ratios of norms of sample points around the circle
						auto tracking_success = CircleTrack(cau_times.front(),cau_samples.front());

						this->IncrementCycleNumber(1);

						if (tracking_success != SuccessCode::Success)
							return tracking_success;

						// find the ratio of the maximum and minimum coordinate wise for the loop. 
						if (RatioEGOperatingZoneTest<CT>())
						{ // then we believe we are in the EG operating zone, since the path is relatively flat.  i still disbelieve this is a good test (dab 20160310)
							while (true)
							{
								if (CheckClosedLoop<CT>())
								{//error is small enough, exit the loop with success. 
									return_value = SuccessCode::Success;
									continue_loop = false;
									break;
								}
								else if(this->CycleNumber() > fail_safe_max_cycle_number)
								{// too many iterations
									return_value = SuccessCode::CycleNumTooHigh;
									continue_loop = false;
									break;
								}

								//compute next loop, the last sample in times and samples is the sample our loop ended on. Either where we started or on another sheet at the same time value. 
								tracking_success = CircleTrack(cau_times.back(),cau_samples.back());

								this->IncrementCycleNumber(1);

								if(tracking_success != SuccessCode::Success)
									return tracking_success;
							}

							if(return_value == SuccessCode::CycleNumTooHigh)
							{//see if we should continue to the next sample point
								if(abs(cau_times.back()) < this->EndgameSettings().min_track_time)
									continue_loop = false;
								else
									continue_loop = true;
							}
						}//end if (RatioEGOperatingZoneTest())
						else 
						{
							//find the time for the next sample point
							next_time = ps_times.back() * RT(this->EndgameSettings().sample_factor);

							SuccessCode tracking_success = this->GetTracker().TrackPath(next_sample,ps_times.back(),next_time,ps_samples.back());
							ps_times.pop_front();
							ps_samples.pop_front();

							ps_times.push_back(next_time);
							ps_samples.push_back(next_sample);

							if(tracking_success != SuccessCode::Success)
								return tracking_success;
						}
					} //end while(continue_loop)

					return return_value;
				}//end InitialCauchyLoops


				/*
				Input: Time and sample at the endgame boundary (usually 0.1). Along with address of time and sample for the approximation computed. 

				Output: SuccessCode deeming if we were successful. 

				Details: This function is in charge of finding the very first approximation of the origin. It does this by first computing some initial samples 
						 like what is done in the Power Series Endgame. We continue to track forward in this manner until we have stabilization of the cycle number being approximated. 
						 This prevents the unnecessary circle tracking if we are possibly not in the endgame operating zone. 

						 Once we have stabilization we then perform InitialCauchyLoops while get the accurate cycle number, and check the norms of the samples and make sure we are ready 
						 to approximate. When ready we call ComputePSEGApproximationOfXAtT0. This function will use a hermtie interpolater to get an approximation of the value at the origin. 
				*/
				/**
				\brief The Cauchy endgame will first find an initial approximation using the notion of the power series endgame. This function computes this approximation and returns a
				SuccessCode to let us know if an error was encountered. 

				\param[out] start_time The time value at which we start the endgame. 
				\param start_point The current space point at start_time.
				\param approximation_time the time at which we want to compute an approximation, usually the origin. 
				\param approximation The space value at approximation_time that we are computing. 

				\tparam CT The complex number type.
				*/

				template<typename CT>
				SuccessCode InitialPowerSeriesApproximation(CT const& start_time, Vec<CT> const& start_point, CT & approximation_time, Vec<CT> & approximation)
				{	
					using RT = typename Eigen::NumTraits<CT>::Real;

					//initialize array holding c_over_k estimates
					std::deque<RT> c_over_k; 

					auto& ps_times = std::get<TimeCont<CT> >(pseg_times_);
					auto& ps_samples = std::get<SampCont<CT> >(pseg_samples_);

					//Compute initial samples for pseg
					auto initial_sample_success = this->ComputeInitialSamples(start_time, start_point, ps_times, ps_samples);
					if (initial_sample_success!=SuccessCode::Success)
						return initial_sample_success;

					c_over_k.push_back(ComputeCOverK<CT>());

					Vec<CT> next_sample;
					CT next_time;

					unsigned ii = 0;
					//track until for more c_over_k estimates or until we reach a cutoff time. 
					while ( (ii < cauchy_settings_.num_needed_for_stabilization) )
					{	
						next_time = ps_times.back() * RT(this->EndgameSettings().sample_factor);

						auto tracking_success = this->GetTracker().TrackPath(next_sample,ps_times.back(),next_time,ps_samples.back());
						if (tracking_success!=SuccessCode::Success)
							return tracking_success;

						ps_samples.pop_front();
						ps_times.pop_front();
						ps_samples.push_back(next_sample);
						ps_times.push_back(next_time);
						c_over_k.push_back(ComputeCOverK<CT>());

						++ii;
					}//end while


					//have we stabilized yet? 
					while(!CheckForCOverKStabilization(c_over_k) && abs(ps_times.back()) > cauchy_settings_.cycle_cutoff_time)
					{
						next_time = ps_times.back() * RT(this->EndgameSettings().sample_factor);

						auto tracking_success = this->GetTracker().TrackPath(next_sample,ps_times.back(),next_time,ps_samples.back());
						if(tracking_success != SuccessCode::Success)
							return tracking_success;

						c_over_k.pop_front();
						ps_samples.pop_front();
						ps_times.pop_front();

						ps_samples.push_back(next_sample);	
						ps_times.push_back(next_time);
						c_over_k.push_back(ComputeCOverK<CT>());

					}//end while
	
					auto cauchy_loop_success = InitialCauchyLoops<CT>();
					if(cauchy_loop_success != SuccessCode::Success)
						return cauchy_loop_success;

					return ComputePSEGApproximationAtT0(approximation, approximation_time);

				}//end InitialPowerSeriesApproximation

				/*
				Input: time_t0 is the time value that we wish to interpolate at.


				Output: A new sample point that was found by interpolation to time_t0.

				Details: This function handles computing an approximation at the origin. 
					We compute the derivatives at the different times and samples. We then make sure all samples are to the same precision before refining them to final tolerance. 
					By InitialCauchyLoops we know what the cycle number is sow we convert derivatives and times to the s-plane where s = t^(1/(cyle number).
					We use the converted times and derivatives along with the samples to do a Hermite interpolation which is found in base_endgame.hpp.
				*/

				/**
				\brief This function takes the pseg_samples_ and pseg_times that have been collected and uses them to compute a Hermite interpolation to the time value time_t0. 


				\param[out] time_t0 is the time value that the power series approximation will interpolate to. 

				\tparam CT The complex number type.
				*/
				template<typename CT>
				SuccessCode ComputePSEGApproximationAtT0(Vec<CT>& result, const CT time_t0)
				{
					using RT = typename Eigen::NumTraits<CT>::Real;

					auto& ps_times = std::get<TimeCont<CT> >(pseg_times_);
					auto& ps_samples = std::get<SampCont<CT> >(pseg_samples_);

					//Compute dx_dt for each sample.
					std::deque< Vec<CT> > pseg_derivatives;
					//initialize the derivative otherwise the first computation in the loop will be wrong. 

					Vec<CT> pseg_derivative = -(this->GetSystem().Jacobian(ps_samples[0],ps_times[0]).inverse())*this->GetSystem().TimeDerivative(ps_samples[0],ps_times[0]);
					for(unsigned ii = 0; ii < this->EndgameSettings().num_sample_points; ++ii)
					{	
						// uses LU look at Eigen documentation on inverse in Eigen/LU.
					 	Vec<CT> pseg_derivative = -(this->GetSystem().Jacobian(ps_samples[ii],ps_times[ii]).inverse())*this->GetSystem().TimeDerivative(ps_samples[ii],ps_times[ii]);
						pseg_derivatives.push_back(pseg_derivative);
					}

					//Checking to make sure all samples are of the same precision.
					unsigned max_precision = 0; 

					for(unsigned ii = 0; ii < this->EndgameSettings().num_sample_points;++ii)
						if(Precision(ps_samples[ii](0)) > max_precision)
							max_precision = Precision(ps_samples[ii](0));

					this->ChangePrecision(max_precision);

					// for(unsigned ii = 0; ii < this->EndgameSettings().num_sample_points;++ii)
					// 	if(Precision(ps_samples[ii](0)) < max_precision)
					// 		for(unsigned jj = 0; jj < this->GetSystem().NumVariables();++jj)
					// 			ps_samples[ii](jj).precision(max_precision);


			 		//Conversion to S-plane.
					std::deque<CT> s_times;
					std::deque< Vec<CT> > s_derivatives;



					for(unsigned ii = 0; ii < ps_samples.size(); ++ii)
					{
						s_derivatives.push_back(pseg_derivatives[ii]*(RT(this->CycleNumber())*pow(ps_times[ii],(RT(this->CycleNumber()) - RT(1))/this->CycleNumber())));
						s_times.push_back(pow(ps_times[ii],RT(1)/RT(this->CycleNumber())));

						auto refine_code = this->RefineSample(ps_samples[ii],ps_samples[ii],ps_times[ii]);
						if (refine_code != SuccessCode:: Success)
							return refine_code;
					}

					result = bertini::tracking::endgame::HermiteInterpolateAndSolve(time_t0, this->EndgameSettings().num_sample_points, s_times, ps_samples, s_derivatives);
					return SuccessCode::Success;
				}//end ComputePSEGApproximationOfXAtT0

				/*
				Input: A deque of the cauchy samples needed to compute an approximation using the Cauchy Integral Formula. 

				Output: The approximation computed by using the samples collected from CircleTrack to sum and divide by the total number used to track back to the starting position. 

				Details: We can compute the Cauchy Integral Formula in this particular instance by computing the mean of the samples we have collected around the origin. 
				*/

				/**
				\brief Function that computes the mean of the samples that were collected while tracking around the origin. This value is the approximation of the value at the origin. 

				\param[out] Deque of cauchy_samples, used to help template the function with CT.

				\tparam CT The complex number type.
				*/
				template<typename CT>
				SuccessCode ComputeCauchyApproximationOfXAtT0(Vec<CT>& result)
				{	using RT = typename Eigen::NumTraits<CT>::Real;
					auto& cau_times = std::get<TimeCont<CT> >(cauchy_times_);
					auto& cau_samples = std::get<SampCont<CT> >(cauchy_samples_);

					if (cau_samples.size() != this->CycleNumber() * this->EndgameSettings().num_sample_points+1)
					{
						std::stringstream err_msg;
						err_msg << "to compute cauchy approximation, cauchy_samples must be of size " << this->CycleNumber() * this->EndgameSettings().num_sample_points+1 << " but is of size " << cau_samples.size() << '\n';
						throw std::runtime_error(err_msg.str());
					}

					result = Vec<CT>::Zero(this->GetSystem().NumVariables());//= (cau_samples[0]+cau_samples.back())/2; 

					for(unsigned int ii = 0; ii < this->CycleNumber() * this->EndgameSettings().num_sample_points; ++ii)
					{
						auto refine_code = this->RefineSample(cau_samples[ii],cau_samples[ii],cau_times[ii]);
						if (refine_code!=SuccessCode::Success)
							return refine_code;
						result += cau_samples[ii];
					}
					result /= RT(this->CycleNumber() * this->EndgameSettings().num_sample_points);
					return SuccessCode::Success;

				}


				/*
				Input: A time and sample value to start at for circle tracking around the origin. 

				Output: A SuccessCode deeming if we were successful or not. 

				Details: This function populates the deque cauchy_samples and cauchy_times. These are data members of the class and are not passed in. This function will continue to 
				call CircleTrack until we have closed the loop. 
				*/

				/**
				\brief Function that will utilize CircleTrack and CheckClosedLoop to collect all samples while tracking around the origin till we close the loop. 


				\param[out] starting_time The time value at which we start finding cauchy samples 
				\param starting_sample The current space point at starting_time, this will also be the first cauchy sample. 

				\tparam CT The complex number type.
				*/
				template<typename CT>
				SuccessCode ComputeCauchySamples(CT const& starting_time, Vec<CT> const& starting_sample)
				{
					auto& cau_times = std::get<TimeCont<CT> >(cauchy_times_);
					auto& cau_samples = std::get<SampCont<CT> >(cauchy_samples_);

					cau_times.clear();
					cau_samples.clear();
					cau_times.push_back(starting_time);
					cau_samples.push_back(starting_sample);
					this->CycleNumber(0);


					while( this->CycleNumber() < cauchy_settings_.fail_safe_maximum_cycle_number )
					{
						//track around the origin once.
						auto tracking_success = CircleTrack(cau_times.back(),cau_samples.back());
						this->IncrementCycleNumber(1);

						if(tracking_success != SuccessCode::Success)
						{
							std::cout << "Cauchy loop fail tracking\n\n";
							return tracking_success;
						}
						else if(CheckClosedLoop<CT>())
						{
							return SuccessCode::Success;
						}
					} 

					return SuccessCode::CycleNumTooHigh;
				}//end ComputeCauchySamples
				

				/*
				Input: A time and sample to start the cauchy endgame.

				Output: SuccessCode deeming if we were successful or not.

				Details: This function runs the entire Cauchy Endgame. We first take our endgame boundary time value and sample to find a first approximation of the origin. This is done by
				using the idea for the power series endgame. We check for stabilization of the cycle number, and check to see when the ratios of the maximum and minimum norm of samples collected
				by CircleTrack are withing a tolerance. When both of these conditions are met we do a Hermite interpolation. 

				At this point we can start tracking in to the origin while using CircleTrack to compute samples and calculating their mean to get an approximation of the origin using the Cauchy
				Integral Formula. 
				*/

				/**
				\brief Function that runs the Cauchy endgame.

				To begin, this function will compute a first approximation using the power series endgame notion. This approximation is made after a heuristic on the stabilization of 
				the cyle number is made, and after the maximum and minimum norms of tracked space values around the origin are withing a certain tolerance. 			

				\param[out] start_time The time value at which we start the endgame. 
				\param start_point The current space point at start_time.

				\tparam CT The complex number type.
				*/
				template<typename CT>
				SuccessCode CauchyEG(CT start_time, Vec<CT> start_point)
				{	
					if (this->preserve_precision_)
						this->initial_precision_ = mpfr_float::default_precision();


					using RT = typename Eigen::NumTraits<CT>::Real;

					auto& cau_times = std::get<TimeCont<CT> >(cauchy_times_);
					auto& cau_samples = std::get<SampCont<CT> >(cauchy_samples_);
					auto& ps_times = std::get<TimeCont<CT> >(pseg_times_);
					auto& ps_samples = std::get<SampCont<CT> >(pseg_samples_);


					ClearTimesAndSamples<CT>(); //clear times and samples before we begin.
					this->CycleNumber(0);
					CT origin(0,0);

					CT next_time;
					Vec<CT> next_sample;

					Vec<CT> prev_approx, latest_approx;
					RT approximate_error;
					Vec<CT>& final_approx = std::get<Vec<CT> >(this->final_approximation_at_origin_);
				

					//Compute the first approximation using the power series approximation technique. 
					auto initial_ps_success = InitialPowerSeriesApproximation(start_time, start_point, origin, prev_approx);  // last argument is output here
					if(initial_ps_success != SuccessCode::Success)
						return initial_ps_success;

					next_time = ps_times.back();
					RT norm_of_dehom_of_prev_approx, norm_of_dehom_of_latest_approx;

					if(this->SecuritySettings().level <= 0)
						norm_of_dehom_of_prev_approx = this->GetSystem().DehomogenizePoint(prev_approx).norm();
					do
					{
						//Compute a cauchy approximation.
						auto extrapolation_success = ComputeCauchyApproximationOfXAtT0<CT>(latest_approx);
						if (extrapolation_success!=SuccessCode::Success)
							return extrapolation_success;

						if(this->SecuritySettings().level <= 0)
							norm_of_dehom_of_latest_approx = this->GetSystem().DehomogenizePoint(latest_approx).norm();

						approximate_error = (latest_approx - prev_approx).norm(); //Calculate the error between approximations. 

						// dehom of prev approx and last approx not used because they are not updated with the most current information. However, prev approx and last approx are 
						// the most current. 

						if(approximate_error < this->Tolerances().final_tolerance)
						{
							final_approx = latest_approx;
							return SuccessCode::Success;
						}
						else if(abs(cau_times.front()) < this->EndgameSettings().min_track_time)
						{//we are too close to t = 0 but we do not have the correct tolerance - so we exit
							final_approx = latest_approx;
							return SuccessCode::FailedToConverge;
						}
						else if(this->SecuritySettings().level <= 0 && 
						   norm_of_dehom_of_prev_approx   > this->SecuritySettings().max_norm &&  
						   norm_of_dehom_of_latest_approx > this->SecuritySettings().max_norm  )
						{//we are too large, break out of loop to return error.
							final_approx = latest_approx;
							return SuccessCode::SecurityMaxNormReached;
						}


						prev_approx = latest_approx;
						norm_of_dehom_of_prev_approx = norm_of_dehom_of_latest_approx;

						next_time *= RT(this->EndgameSettings().sample_factor);

						auto tracking_success = this->GetTracker().TrackPath(next_sample,cau_times.back(),next_time,cau_samples.front());
						if (tracking_success != SuccessCode::Success)
							return tracking_success;

						ps_times.push_back(next_time);  ps_times.pop_front();
						ps_samples.push_back(next_sample); ps_samples.pop_front();

						auto cauchy_samples_success = ComputeCauchySamples(next_time,next_sample);

						//Added because of Griewank osborne test case where tracker returns GoingToInfinity.  
						//The cauchy endgame should stop instead of attempting to continue. 
						//This is because the tracker will not track to any meaningful space values. 
						if (cauchy_samples_success == SuccessCode::GoingToInfinity)
							return SuccessCode::GoingToInfinity;
						else if (cauchy_samples_success != SuccessCode::Success && abs(cau_times.front()) < this->EndgameSettings().min_track_time)
						{// we are too close to t = 0 but we do have the correct tolerance -so we exit.
							final_approx = latest_approx;
							return SuccessCode::MinTrackTimeReached;
						}
						else if(cauchy_samples_success != SuccessCode::Success)
						{
							final_approx = latest_approx;
							return cauchy_samples_success;
						}

					} while (approximate_error > this->Tolerances().final_tolerance);

					final_approx = latest_approx;
					return SuccessCode::Success;
				} //end CauchyEG
			};


			
		}//namespace endgame

	}//namespace tracking

} // namespace bertini
#endif
