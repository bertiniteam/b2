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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire
//

/**
\file pool.hpp 
*/


#ifndef BERTINI_GENERIC_POOL_HPP
#define BERTINI_GENERIC_POOL_HPP

namespace bertini {

	namespace detail {
	

	/// \brief Pool pointer policy that holds objects via std::shared_ptr.
	template<typename ObjT>
	struct SharedPointerPolicy
	{
		using PointerType = std::shared_ptr<ObjT>;  ///< The held pointer type.

		/// \brief Make a held pointer from an object value.
		static
		PointerType FromObj(ObjT d)
		{
			return std::make_shared<ObjT>(std::move(d));
		}

		/// \brief Make a held pointer to a default-constructed object.
		static
		PointerType DefaultConstructed()
		{
			return std::make_shared<ObjT>();
		}

	};


	/// \brief The default pool pointer policy (shared_ptr).
	template<typename ObjT>
	using DefaultPointerPolicy = SharedPointerPolicy<ObjT>;



	/// \brief A process-wide cache of objects, held by pointer, that hands out shared references.
	template<typename ObjT, class PointerPolicy = DefaultPointerPolicy<ObjT> >
	class Pool
	{




	public:

		using HeldType = typename PointerPolicy::PointerType;  ///< The pointer type the pool holds.
		using PoolHolderType = std::vector< HeldType >;  ///< The container of held pointers.

		/// \brief Add an object (by value) to the pool, returning the held pointer.
		static
		HeldType NonPtrAdd(ObjT d)
		{
			// held_data_.push_back(PointerPolicy::FromObj(d));
			held_data_.push_back(PointerPolicy::FromObj(d));
			return held_data_.back();
		}

		/// \brief Add an already-held pointer to the pool, returning it.
		static
		HeldType PtrAdd(HeldType d)
		{
			held_data_.push_back(d);
			return held_data_.back();
		}

		/// \brief Construct an object in place in the pool, returning the held pointer.
		template<typename ... Ts>
		static
		HeldType Make(Ts&& ...ts )
		{
			auto thing = std::shared_ptr<ObjT>(new ObjT(ts...));
			held_data_.push_back(thing);
			return held_data_.back();
		}

		/// \brief Drop pooled objects that are no longer referenced anywhere else.
		static
		void PurgeCache()
		{
			std::remove_if(held_data_.begin(), held_data_.end(), [](HeldType const& h){return h.use_count()==1;});
		}

	private:
		static
		PoolHolderType held_data_;  ///< The pool's held pointers.
	};

	/// \cond POOL_STATIC_DEFINITION
	template<typename ObjT, class Policy> typename Pool<ObjT,Policy>::PoolHolderType Pool<ObjT,Policy>::held_data_;
	/// \endcond
	} // re: detail
} // re: bertini

#endif
