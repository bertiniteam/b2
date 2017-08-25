//This file is part of Bertini 2.
//
//pool.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//pool.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with pool.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Dani Brake
// University of Notre Dame
//

/**
\file pool.hpp 
*/


#ifndef BERTINI_GENERIC_POOL_HPP
#define BERTINI_GENERIC_POOL_HPP

namespace bertini {

	namespace detail {
	

	template<typename T>
	struct SharedPointerPolicy
	{
		using PointerType = std::shared_ptr<T>;

		static
		PointerType FromObj(T d)
		{
			return std::make_shared<T>(std::move(d));
		}

		static
		PointerType DefaultConstructed()
		{
			return std::make_shared<T>();
		}

	};


	template<typename T>
	using DefaultPointerPolicy = SharedPointerPolicy<T>;



	template<typename T, class PointerPolicy = DefaultPointerPolicy<T> >
	class Pool
	{
		using HeldType = typename PointerPolicy::PointerType;

		using PoolHolderType = std::vector< HeldType >;

		PoolHolderType held_data_;

	public:


		HeldType NonPtrAdd(T d)
		{
			held_data_.push_back(PointerPolicy::FromObj(d));
			return held_data_.back();
		}

		HeldType PtrAdd(HeldType d)
		{
			held_data_.push_back(d);
			return held_data_.back();
		}

		HeldType NewObj()
		{
			held_data_.push_back(PointerPolicy::DefaultConstructed());
			return held_data_.back();
		}

		void PurgeCache()
		{
			std::remove_if(held_data_.begin(), held_data_.end(), [](HeldType const& h){return h.use_count()==1;});
		}

	};

	} // re: detail
} // re: bertini

#endif
