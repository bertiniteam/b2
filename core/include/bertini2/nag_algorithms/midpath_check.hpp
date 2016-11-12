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
			using PathIndT = unsigned long long;
			
			Midpath( StartSystemT start_system, RealType near_tolerance)
			{
				boundary_near_tolerance_ = near_tolerance;
				
				start_system_ = start_system;
				
				passed_ = true;
			}
			
			
			/**
			 \struct Stores data for each path that crosses, including the index of the path, all paths that it crosses with, and whether it has the same
				starting point as a path it crosses with.
			 */
			struct CrossedPath{
				
				CrossedPath(PathIndT index, PathIndT crossed_with, bool same_start)
				{
					index_ = index;
					crossed_with_.push_back( std::make_pair(crossed_with, same_start) );
					rerun_ = !same_start;
				}
				
				
				PathIndT index() const
				{
					return index_;
				}
				
				void crossed_with(std::pair<PathIndT, bool> crossed)
				{
					crossed_with_.push_back(crossed);
				}
				
				void rerun(bool same_start)
				{
					rerun_ = rerun_ || !same_start;
				}
				
				bool rerun() const
				{
					return rerun_;
				}


			private:
				PathIndT index_;
				std::vector< std::pair<PathIndT, bool> > crossed_with_; // first is index of path crossed with, second is bool = true if had same starting point as crossed path
				bool rerun_;


			}; //re: struct CrossedPath


			
			bool Passed()
			{
				return passed_;
			}
			
			
			
			
			/**
			 \brief Checks the solution data at the endgame boundary to see if any paths have crossed during tracking before the endgame.
			 
			 \param boundary_data Solution data at the endgame boundary
			 
			 \returns A Data object which stores data about crossed paths
			 
			*/
			bool Check(BoundaryData const& boundary_data)
			{
				for (PathIndT ii = 0; ii < boundary_data.size(); ++ii)
				{
					const Vec<ComplexType>& solution_ii = std::get<Vec<ComplexType>>(boundary_data[ii]);
										
					for(PathIndT jj = ii+1; jj < boundary_data.size(); ++jj)
					{
						if(std::get<tracking::SuccessCode>(boundary_data[ii]) == tracking::SuccessCode::Success)
						{
							const Vec<ComplexType>& solution_jj = std::get<Vec<ComplexType>>(boundary_data[jj]);
							const Vec<ComplexType> diff_sol = solution_ii - solution_jj;
							
							if((diff_sol.template lpNorm<Eigen::Infinity>()/solution_ii.template lpNorm<Eigen::Infinity>()) < boundary_near_tolerance_)
							{
								bool i_already_stored = false;
								bool j_already_stored = false;
								// Check if start points are the same
								const auto& start_ii = start_system_.template StartPoint<ComplexType>(ii);
								const auto& start_jj = start_system_.template StartPoint<ComplexType>(jj);
								auto diff_start = start_ii - start_jj;
								bool same_start = (diff_start.norm() > boundary_near_tolerance_);
								
								
								// Check if path has already been stored in crossed_paths_
								for(auto& v : crossed_paths_)
								{
									if(v.index() == ii)
									{
										v.crossed_with(std::make_pair(jj,same_start));
										i_already_stored = true;
										v.rerun(same_start);
									}
									if(v.index() == jj)
									{
										v.crossed_with(std::make_pair(ii,same_start));
										j_already_stored = true;
										v.rerun(same_start);
									}
								}
								
								// If not already stored, create CrossedPath object and add to crossed_paths_
								if(!i_already_stored)
								{
									CrossedPath tempPath(ii, jj, same_start);
									crossed_paths_.push_back(tempPath);
								}
								if(!j_already_stored)
								{
									CrossedPath tempPath(jj, ii, same_start);
									crossed_paths_.push_back(tempPath);
								}
								
								passed_ = false;
								
							}
						}
					}
				}
				return passed_;
			};
			
			
			std::vector<CrossedPath> GetCrossedPaths()
			{
				return crossed_paths_;
			}
			
			
			
		private:
			RealType boundary_near_tolerance_;
			StartSystemT start_system_;
			
			std::vector<CrossedPath> crossed_paths_; // Data for all paths that crossed on last check
			bool passed_;  // Did the check pass?
			
			
		}; //re: struct Midpath

	} // re: namespace algorithm
}// re: namespace bertini
