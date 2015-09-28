//This file is part of Bertini 2.0.
//
//powerseries.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//powerseries.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with predict.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  powerseries.hpp
//
//  copyright 2015
//  Tim Hodges
//  Colorado State University
//  Department of Mathematics
//  Fall 2015
#ifndef BERTINI_TRACKING_POWERSERIES_HPP
#define BERTINI_TRACKING_POWERSERIES_HPP
#include <queue>
#include "tracking/step.hpp"

#include "system.hpp"
namespace bertini{ 

	namespace tracking {

		namespace endgame {
		/**
		This is a header for all functions need to perform the powerseries endgame. 

		\param target_time is the time for which we want to approximate using our interpolation.
		\param Samples a queue of points in space that will be used in the interpolation.
		\param Derivatives is a queue of derivatives in the same order of the points in Samples. 
				The order of these derivatives are in order of when they are computed.
				dH_dx, dH_dp, dp_dt are readily available. 
				from the above we get dH_dt, via. dH_dt = dH_dx * dx_dt + dH_dp * dp_dt = 0
				Then from dH_dt we obtain dx_dt = -(dH_dx)^(-1) * (dH_dp * dp_dt)
				t = s^c --> s = t^(1/c) 
				dt_ds = c*s^(c-1) = c*(t)^((c-1)/c)
				lastly, dx_ds = dx_dt * dt_ds

				So the ordering is [dH_dx, dH_dp, dp_dt, dH_dt, dx_dt, dt_ds, dx_ds]



		\param sys The system being solved.
		\param current_space The current space variable vector.
		\param current_time The current time.
		\param delta_t The size of the time step.
		\param condition_number_estimate The computed estimate of the condition number of the Jacobian.
		\param num_steps_since_last_condition_number_computation.  Updated in this function.
		\param frequency_of_CN_estimation How many steps to take between condition number estimates.
		\param prec_type The operating precision type.  
		\param tracking_tolerance How tightly to track the path.
		\param AMP_config The settings for adaptive multiple precision.

		\tparam ComplexType The complex number type for evaluation.
		\tparam RealType The complex number type for evaluation.
		*/



			template<typename ComplexType>
			Vec<ComplexType> HermiteInterpolateAndSolve(ComplexType const& target_time, 
				std::queue< Vec<ComplexType> > const& samples, 
				std::queue< Vec<ComplexType> > const& derivatives){
				//Hermite Interpolation needs to be implemented	


			}

			template<typename ComplexType, typename RealType>
			unsigned int BoundOnCycleNumber(ComplexType sample2, ComplexType sample1, ComplexType sample0,
				RealType sample_factor){
				Vec<ComplexType> rand_vector = Vec<ComplexType>::Random(sample0.size);

				// need to write the code for estimate. Make sure to use Eigen and transpose to use LinAlg
				// need to have cycle_number_amplification, this is the 5 used for 5 * estimate


			}

			template<typename ComplexType>
			unsigned int ComputeCycleNumber(ComplexType current_time, Vec<ComplexType> x_current_time,
				unsigned int upper_bound_on_cycle_number,
				std::queue< Vec<ComplexType> >const& samples,
				std::queue< Vec<ComplexType> >const& derivatives){
				/*
					min_diff_in_approximations = 1e300;
					min_found_difference = min_diff_in_approximations
					cycle_number = 0;

					for c = 1:upper_bound_on_cycle_number
						s_current_time = pow(current_time,1/c);
						dx_ds = dx_dt * c *(pow(t,((c-1)/c));
						approx = hermite....
						if norm < min_found_difference 
							min_found_difference = norm
							cycle_number = c;

						end

					end 

					return cycle_number 

				*/	

			}


			template<typename ComplexType, typename RealType>
			std::queue< Vec<ComplexType> > ComputeInitialSamples(ComplexType endgame_time,
				Vec<ComplexType> x_endgame, 
				unsigned int num_samples,
				RealType sample_factor){

				/*
					This will find the first num_samples samples we need to make power series endgame run. 
					First sample is (x_endgame, t_endgame)
					push first sample on to samples queue. 

					for i = 2:num_samples
						next_time = current_time * sample_factor;
						new_sample = track(current_time, x_current_time, next_time)
						samples.push_back <new_sample, next_time>
						current_time = next_time
						x(current_time) = new_sample
						

					end
					return samples
				*/	


			}


			template<typename ComplexType, typename RealType>
			Vec<ComplexType> ComputeApproximationOfXAtT0(ComplexType time_t0, 
				std::queue< Vec<ComplexType> > const& samples,
				System & sys){

				//Steps: check all samples are at the same precision first	
				//	2. calculate dx of H @ samples
				//  upper_bound = BoundOnCycleNumber
				//	cycle_number = ComputeCycleNumber 
				//  Conversion to S Plane
				//  appproximation_of_x_at_t0 = HermiteInterpolateAndSolve

			}

			template<typename ComplexType, typename RealType>
			Vec<ComplexType> PSEG(ComplexType endgame_time, Vec<ComplexType> x_current_time,
			 System & sys){

				//using num_samples = bertini::num_samples;
				unsigned int num_samples;
				//using sample_factor = bertini::sample_factor;
				RealType sample_factor;
				//using final_tolerance = bertini::final_tolerance;
				RealType final_tolerance;
				//using min_track_time = bertini::min_track_time;
				RealType min_track_time;
				//using mpfr_float = boost::multiprecision::mpfr_float;



			 	std::queue< Vec<ComplexType> > samples = ComputeInitialSamples(endgame_time, x_current_time, num_samples,
			 		sample_factor);

			 	RealType norm_last = 0;

			 	

			 	ComplexType	origin = mpfr_float("0","0");
			 	

			 	Vec<ComplexType> prev_approx = ComputeApproximationOfXAtT0(origin, samples, sys);

			 	prev_approx = dehom(prev_approx);
			 	
			 	norm_last = norm(prev_approx);

			 	RealType approx_error = 1;


			 	while (approx_error < final_tolerance){

			 		ComplexType current_time = current_time * sample_factor;

			 		if (norm(current_time) < min_track_time){
			 			//error as you are past the minimum you are allowed to track to.
			 			//return
			 		}

			 		Vec<ComplexType> x_current_time = track(current_time, current_time/sample_factor, x_current_time);

			 		// if(norm_last > bertini::SecurityMaxNorm){
			 		// 	//error and return
			 		// }

			 		samples.push(x_current_time);
			 		samples.pop();

			 		//push x_current_time into queue

			 		Vec<ComplexType> latest_approx = ComputeApproximationOfXAtT0(endgame_time,samples, sys);
			 		//latest_approx = ComputeApproximationOfAtT0

			 		approx_error = norm(latest_approx - prev_approx);
				//	approx_error = norm(latest_approx, prev_approx)

			 		if(approx_error < final_tolerance){
			 			printf("YES!!");
			 			return latest_approx;
			 		}

					// if (approx_error < final_tolerance){
					// 	SUCCESS 
					// 	return
					// 	}

					prev_approx = latest_approx;

					norm_last = inf_norm(dheom(prev_approx));
							

			 	} //end while	

				
					

			} //end PSEG

		}//namespace endgame

	}//namespace tracking

} // namespace bertini
#endif
