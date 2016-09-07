//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/midpath_check.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/midpath_check.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/midpath_check.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// James Collins, West Texas A&M University


/**
\file bertini2/nag_algorithms/midpath_check.hpp 

\brief Provides methods for checking for path crossings.  

This essentially amounts to being certain that no two points are the same.
*/


#pragma once

#include"bertini2/nag_algorithms/config.hpp"
#include "bertini2/tracking.hpp"


namespace bertini{
	namespace algorithm{

		template<typename StartSystemT, typename RealType, typename ComplexType>
		struct Midpath
		{
			
			using BoundaryData = SolnCont< std::tuple<Vec<ComplexType>, tracking::SuccessCode, RealType>>;
			
			Midpath( StartSystemT start_system, RealType near_tolerance)
			{
				boundary_near_tolerance_ = near_tolerance;
				
				start_system_ = start_system;
			}
			
			struct Data{
				bool Passed()
				{
					return pass_;
				}
				
				void Crossed(unsigned long long path_index, bool same_start)
				{
					bool already_stored = false;
					for(auto v : crossed_paths)
					{
						auto index = std::get<unsigned long long>(v);
						if(path_index == index)
						{
							already_stored = true;
						}
					}
					
					if(already_stored == false)
					{
						crossed_paths.push_back(std::make_tuple(path_index, same_start));
					}
					
					pass_ = false;
				}
				
				auto GetCrossedPaths() const
				{
					return crossed_paths;
				}


			private:
				bool pass_ = true;
				std::vector< std::tuple< unsigned long long, bool > > crossed_paths;


			}; //re: struct Data


			
			
			
			
			
			/**
			 \brief Checks the solution data at the endgame boundary to see if any paths have crossed during tracking before the endgame.
			 
			 \param boundary_data Solution data at the endgame boundary
			 
			 \returns A Data object which stores data about crossed paths
			 
			*/
			Data Check(BoundaryData const& boundary_data)
			{
				Data check_data;
				for (unsigned long long ii = 0; ii < boundary_data.size(); ++ii)
				{
					const Vec<ComplexType>& solution_ii = std::get<Vec<ComplexType>>(boundary_data[ii]);
										
					for(unsigned long long jj = ii+1; jj < boundary_data.size(); ++jj)
					{
						if(std::get<tracking::SuccessCode>(boundary_data[ii]) == tracking::SuccessCode::Success)
						{
							const Vec<ComplexType>& solution_jj = std::get<Vec<ComplexType>>(boundary_data[jj]);
							const Vec<ComplexType> diff_sol = solution_ii - solution_jj;
							
							if(diff_sol.norm() < boundary_near_tolerance_)
							{
								// Check if start points are the same
								const auto& start_ii = start_system_.template StartPoint<ComplexType>(ii);
								const auto& start_jj = start_system_.template StartPoint<ComplexType>(jj);
								auto diff_start = start_ii - start_jj;
								bool same_start = (diff_start.norm() > boundary_near_tolerance_);
								check_data.Crossed(ii, same_start);
								check_data.Crossed(jj, same_start);
							}
						}
					}
				}
				return check_data;
			};
			
			
		private:
			RealType boundary_near_tolerance_;
			StartSystemT start_system_;
		}; //re: struct Midpath

	} // re: namespace algorithm
}// re: namespace bertini
