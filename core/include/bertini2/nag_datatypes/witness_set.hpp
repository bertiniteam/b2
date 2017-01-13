//This file is part of Bertini 2.
//
//bertini2/nag_datatypes/witness_set.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_datatypes/witness_set.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_datatypes/witness_set.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/nag_datatypes/witness_set.hpp 

\brief Provides the witness set data type for Bertini2.
*/


#pragma once

#include <deque>
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/nag_datatypes/common/policies.hpp"
#include "bertini2/system/slice.hpp"
#include "bertini2/system/system.hpp"


namespace bertini {

	namespace nag_datatype {

		template <typename T> 
		using PointCont = std::deque<T>;

		using IndexT = unsigned long long;

		/**
		\brief The main way to represent a positive dimensional component in numerical algebraic geometry.
	

		A witness set consists of three properties:
	
		1. Points
		2. A generic linear slice
		3. The system.  

		The witness set captures two of the most important properties to know about a component of an algebraic variety: 

		 * degree
		 * dimension

		 The degree is represented in terms of points.  The number of points in the witness set is the degree of the component.

		 The dimension is represented in terms of either the slice, or the system itself.  The number of individual linears in the total slice should be equal to the difference of the number of variables in the system, minus the number of functions.  That is, the system's underdeterminedness is the dimension.  This assumes that randomization has been done on incomplete intersections -- precisely those systems for which this equality fails.

		 \see Randomize
		*/
		template<typename NumT, typename SystemT = System, template<typename> class ObjManagementP = policy::CopyGiven >
		class WitnessSet
		{
			using PointP = ObjManagementP<Vec<NumT>>;
			using SliceP = ObjManagementP<LinearSlice>;
			using SystemP = ObjManagementP<SystemT>;

			using PointContT = PointCont<typename PointP::HeldT>;

			PointContT points_;
			typename SliceP::HeldT slice_;
			typename SystemP::HeldT system_;


public:

			WitnessSet(){}

			WitnessSet(PointContT const& pts, LinearSlice const& slc, System const& sys) :
				points_(pts),
				slice_(SliceP::AtSet(slc)),
				system_(SystemP::AtSet(sys))
			{}

			void SetSlice(LinearSlice const& s)
			{
				slice_ = SliceP::AtSet(s);
			}

			const LinearSlice & GetSlice() const
			{
				return SliceP::AtGet(slice_);
			}

			inline
			void AddPoint(Vec<NumT> const& p)
			{
				points_.push_back(PointP::AtSet(p));
			}

			void SetPoints(PointContT && pts)
			{
				points_ = PointP::AtSet(pts);
			}

			inline
			const PointContT& GetPoints() const
			{
				return points_;
			}


			const Vec<NumT>& GetPoint(IndexT ind) const
			{
				if (ind > Degree())
				{
					std::stringstream ss;
					ss << "asking for point of index " << ind << ", outside range [0,"<<Degree()-1 << "].";
					throw std::runtime_error(ss.str());
				}

				return points_[ind];
			}

			inline
			auto Degree() const
			{
				return points_.size();
			}

			inline
			auto Dimension() const
			{
				return slice_.Dimension();
			}

			void SetSystem(SystemT && sys)
			{
				system_ = SystemP::AtSet(sys);
			}

			inline
			const SystemT& GetSystem() const
			{
				return system_;
			}


			/**
			\brief Check whether the witness set makes sense, in terms of system size, number of variables, and linear slice size.
			*/
			bool IsConsistent() const
			{
				return (system_.NumVariables() - system_.NumFunctions()) == (slice_.Dimension());
			}


		};
	}
}


