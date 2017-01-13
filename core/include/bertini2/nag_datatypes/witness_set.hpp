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


		template<typename NumT, template<typename> class ObjManagementP = policy::CopyGiven >
		class WitnessSet
		{
			using PointP = ObjManagementP<std::deque<Vec<NumT>>>;
			using SliceP = ObjManagementP<LinearSlice>;
			using SystemP = ObjManagementP<System>;


			typename PointP::HeldT points_;
			typename SliceP::HeldT linear_slice_;
			typename SystemP::HeldT system_;


public:

			WitnessSet(){}

			WitnessSet(std::deque<Vec<NumT>> && pts, LinearSlice && slc, System && sys) :
				points_(PointP::AtSet(pts)),
				linear_slice_(SliceP::AtSet(slc)),
				system_(SystemP::AtSet(sys))
			{}

			void SetSlice(LinearSlice const& s)
			{
				linear_slice_ = ObjManagementP<LinearSlice>::AtSet(s);
			}

			const LinearSlice & GetSlice() const
			{
				return ObjManagementP<LinearSlice>::AtGet(linear_slice_);
			}


			void AddPoint(Vec<NumT> const& p)
			{
				points_.push_back(p);
			}


			auto Degree() const
			{
				return points_.size();
			}

			auto Dimension() const
			{
				return linear_slice_.Dimension();
			}
		};
	}
}


