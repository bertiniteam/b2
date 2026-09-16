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
// Copyright(C) Bertini2 Development Team
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
#include <boost/serialization/deque.hpp>
#include <boost/serialization/split_member.hpp>
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/nag_datatypes/common/policies.hpp"
#include "bertini2/system/slice.hpp"
#include "bertini2/system/system.hpp"


namespace bertini {

	namespace nag_datatype {

		/// \brief The container type used to hold points (a deque).
		template <typename T>
		using PointCont = std::deque<T>;

		using IndexT = typename PointCont<int>::size_type;  ///< The index type for point containers.

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
	
		\tparam ObjManagementP  The policy class for object management.  Default is to copy.  Other provided options are reference wrappers, and shared pointers.

		\see bertini::nag_datatype::policy::Copy, bertini::nag_datatype::policy::Reference, bertini::nag_datatype::policy::SharedPtr

		 \see Randomize
		*/
		template<typename NumT, typename SystemT = System, template<typename> class ObjManagementP = policy::Copy >
		class WitnessSet
		{
			using PointP = ObjManagementP<Vec<NumT>>;
			using SliceP = ObjManagementP<Slice>;
			using SystemP = ObjManagementP<SystemT>;

			using PointContT = PointCont<typename PointP::HeldT>;

			PointContT points_;
			typename SliceP::HeldT slice_;
			typename SystemP::HeldT system_;


public:
			/**
			\brief The default constructor.  

			The instantiability of this constructor depends on the default-constructability of the policy class used for held objects.
			*/
			WitnessSet(){}

			/**
			\brief Construct a witness set from points, slice, and system.

			Uses the held type of the policy class for object management.
			*/
			WitnessSet(PointContT const& pts, typename SliceP::HeldT const& slc, typename SystemP::HeldT const& sys) :
				points_(pts),
				slice_(slc),
				system_(sys)
			{}


			/**
			\brief Get the degree of the witness set.  

			This is the number of points in the set.
			*/
			inline
			auto Degree() const
			{
				return GetPoints().size();
			}


			/**
			\brief Query the dimension of the witness set.  Fundamentally, this is the dimension of the slice.
			*/
			inline
			auto Dimension() const
			{
				return GetSlice().Dimension();
			}


			/**
			\brief Adds a point to the witness set.
			*/
			inline
			void AddPoint(typename PointP::HeldT const& p)
			{
				points_.push_back(p);
			}


			/**
			\brief Sets the set of points for the witness set.
			*/
			void SetPoints(PointContT const& pts)
			{
				points_ = pts;
			}


			/**
			\brief Get (a const reference to) the set of points for the witness set.
			*/
			inline
			const PointContT& GetPoints() const
			{
				return points_;
			}

			/**
			\brief Range-checked version of point getter for witness set.
			*/
			const Vec<NumT>& GetPoint(IndexT ind) const
			{
				if (ind > Degree())
				{
					std::stringstream ss;
					ss << "asking for point of index " << ind << ", outside range [0,"<<Degree()-1 << "].";
					throw std::runtime_error(ss.str());
				}

				return PointP::AtGet(points_[ind]);
			}


			/**
			\brief Non-range-checked accessor to get points from witness set.
			*/
			const Vec<NumT>& operator[](IndexT ind) const
			{
				return PointP::AtGet(points_[ind]);
			}


			


			/**
			Sets the system for the witness set.  Uses policy for held type.

			\see bertini::nag_datatype::policy namespace
			*/
			void SetSystem(typename SystemP::HeldT const& sys)
			{
				system_ = sys;
			}


			/**
			\brief Gets a (const reference to) the system for the witness set.
			*/
			inline
			const SystemT& GetSystem() const
			{
				return SystemP::AtGet(system_);
			}


			/**
			\brief Sets the slice for the witness set.

			uses the held type policy.
			*/
			void SetSlice(typename SliceP::HeldT const& s)
			{
				slice_ = s;
			}

			/**
			Gets (a const reference to) the slice for the witness set.
			*/
			const Slice & GetSlice() const
			{
				return SliceP::AtGet(slice_);
			}

			/**
			\brief Check whether the witness set makes sense, in terms of system size, number of variables, and linear slice size.
			*/
			bool IsConsistent() const
			{
				return (GetSystem().NumVariables() - GetSystem().NumNaturalFunctions()) == GetSlice().Dimension();
			}


		private:

			friend class boost::serialization::access;

			/// Serialize the witness set's three parts (points, slice, system).  Save/load are split
			/// so load can restore the system's differentiated form + working precision afterwards:
			/// System::serialize persists neither, so without this a deserialized witness set's system
			/// would be unusable for evaluation/tracking.  Doing it here (not just in the Python pickle
			/// suite) makes the C++ deserialize path -- the one MPI and threading use -- correct too.
			template <typename Archive>
			void save(Archive& ar, const unsigned /*version*/) const
			{
				ar & points_;
				ar & slice_;
				ar & system_;
			}

			template <typename Archive>
			void load(Archive& ar, const unsigned /*version*/)
			{
				ar & points_;
				ar & slice_;
				ar & system_;
				// GetSystem() is const and returns a const System&; Differentiate()/precision() are
				// const (they mutate only the system's mutable evaluation state), so this restores
				// the held system in place regardless of object-management policy.
				GetSystem().Differentiate();
				// the precision self-assignment that used to sit here (a re-materialize at the
				// value already held) is unnecessary: the first evaluation after load
				// materializes the system at the precision of the point it is given (#377)
			}

			BOOST_SERIALIZATION_SPLIT_MEMBER()

		};
	}
}


