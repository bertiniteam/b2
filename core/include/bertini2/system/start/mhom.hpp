//This file is part of Bertini 2.
//
//mhom.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mhom.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with mhom.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Tim Hodges, Colorado State University
// silviana amethyst, university of wisconsin eau claire

/**
\file mhom.hpp 

\brief Defines the MHomogeneous start system type.
*/

#pragma once

#include "bertini2/system/start_base.hpp"
#include "bertini2/system/start/utility.hpp"

namespace bertini 
{
	namespace start_system
	{


		/**
		\brief The $m$-homogeneous start system for Numerical Algebraic Geometry

		
		*/
		class MHomogeneous : public StartSystem
		{
		public:
			MHomogeneous() = default;
			virtual ~MHomogeneous() = default;
			
			/**
			 Constructor for making a multi-homogeneous start system from a polynomial system

			 \throws std::runtime_error, if the input target system is not square, is not polynomial, has a path variable already.

			 Constructor can take in a non-homogeneous system or homogeneous system. 
			*/
			MHomogeneous(System const& s);

			/**
			\brief Creates a degree matrix for constructing the multi-homogeneous start system.
			*/
			void CreateDegreeMatrix(System const& s);

			/**
			\brief Creates all valid partitions for multi-homogeneous start system to create start points.
			*/
			void GenerateValidPartitions(System const& s);

			/**
			\brief Helper function that is used to find valid partitions in the degree matrix.
			*/
			int ChooseColumnInRow(System const& s,Vec<int>& variable_group_counter, int row, int column);

			/**
			Get the number of start points for this m-homogeneous start system.  This is the Bezout bound for the target system.  Provided here for your convenience.
			*/
			unsigned long long NumStartPoints() const override;

			unsigned long long NumStartPointsForPartition(Vec<int> partition) const;

			MHomogeneous& operator*=(Nd const& n);

			MHomogeneous& operator+=(System const& sys) = delete;

			/**
			 \brief Degree matrix holding degrees for all functions in terms of all variable groups
			 */
			Mat<int> degree_matrix_; // stores degrees of all functions in all homogeneous variable groups.
			
			/**
			 \brief Partitions used for creating start points in the multi-homogeneous start system.
			 */
			std::deque< Vec<int> > valid_partitions_;


		private:

			/**
			Get the ith start point, in double precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<complex_dbl> GenerateStartPoint(complex_dbl,unsigned long long index) const override;

			/**
			Get the ith start point, in current default precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<complex_mp> GenerateStartPoint(complex_mp,unsigned long long index) const override;
			
			/**
			 A local version of GenerateStartPoint that can be templated
			*/
			template<typename T>
			void GenerateStartPointT(Vec<T>& start_point, unsigned long long index) const;
			
			

			
			std::vector<unsigned long long> degrees_; ///< stores the degrees of the functions.
			std::vector< VariableGroup > var_groups_;
			/// Random linear-factor coefficients, one matrix per (function, variable group).
			/// Entry (i,j) is a (degree_matrix_(i,j)) x (group_j_size + 1) matrix: each row is
			/// a linear factor over group j's variables, the trailing column being the factor's
			/// constant (0 for projective groups -- their factors are homogeneous).  This is the
			/// data the retired node::LinearProduct used to hold; the products-of-linears block
			/// and the start-point solve read it directly.
			Mat<Mat<complex_mp>> linear_coeffs_;
			std::vector< std::vector<size_t> > variable_cols_; ///< The columns associated with each variable.  The first index is the variable group, the second index is the particular variable in the group.
			size_t num_hom_groups_ = 0; ///< how many of the leading entries of var_groups_ are projective (homogeneous) groups; the rest are affine.  A projective group of size k spans P^{k-1}: dimension k-1, so it takes k-1 functions and its start-point component is a homogeneous vector.

			mutable Vec<complex_mp> temp_v_mp_;

			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned /*version*/) 
			{
				ar & boost::serialization::base_object<StartSystem>(*this);
				ar & degrees_;
			}

		};
	}//end start_system namespace
}//end bertini namespace



