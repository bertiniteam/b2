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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file mhom.hpp 

\brief Defines the MHomogeneous start system type.
*/

#pragma once

#include "bertini2/system/start_base.hpp"
#include "bertini2/limbo.hpp"


namespace bertini 
{
	namespace start_system{


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

			 \throws std::runtime_error, if the input target system is not square, is not polynomial, has a path variable already, or has any homogeneous variable groups.
			*/
			MHomogeneous(System const& s);


			

			/**
			Get the number of start points for this total degree start system.  This is the Bezout bound for the target system.  Provided here for your convenience.
			*/
			unsigned long long NumStartPoints() const override;

			MHomogeneous& operator*=(Nd const& n);

			MHomogeneous& operator+=(System const& sys) = delete;
			
		private:

			/**
			Get the ith start point, in double precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<dbl> GenerateStartPoint(dbl,unsigned long long index) const override;

			/**
			Get the ith start point, in current default precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<mpfr> GenerateStartPoint(mpfr,unsigned long long index) const override;

			std::vector<std::shared_ptr<node::Rational> > random_values_; ///< stores the random values for the start functions.  x^d-r, where r is stored in this vector.
			std::vector<unsigned long long> degrees_; ///< stores the degrees of the functions.


			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<StartSystem>(*this);
				ar & random_values_;
				ar & degrees_;
			}

		};
	}
}



