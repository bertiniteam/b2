//This file is part of Bertini 2.
//
//total_degree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//total_degree.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with total_degree.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire

/**
\file total_degree.hpp 

\brief Defines the TotalDegree start system type.
*/

#pragma once

#include "bertini2/system/start_base.hpp"
#include "bertini2/limbo.hpp"


namespace bertini 
{
	namespace start_system{


		/**
		\brief Total degree start system for 1-homogeneous polynomial systems.

		The most basic and easy-to-construct start system in Numerical Algebraic Geometry.  

		The total degree start system uses functions of the form \f$x_i^{d_i} - r_i\f$, where \f$i\f$ is the index of the function relative to the system, \f$d_i\f$ is the degree of that function, and \f$r_i\f$ is a random complex number.  This is very similar to the roots of unity, except that they are moved away from being centered around the origin, to being centered around a random complex number.  

		Note that the corresponding target system MUST be square -- have the same number of functions and variables.  The start system cannot be constructed otherwise, particularly because it is written to throw at the moment if not square.

		The start points are accesses by index (unsigned long long), instead of being generated all at once.
		*/
		class TotalDegree : public StartSystem
		{
		public:
			TotalDegree() = default;
			virtual ~TotalDegree() = default;
			
			/**
			 Constructor for making a total degree start system from a polynomial system

			 \throws std::runtime_error, if the input target system is not square, is not polynomial, has a path variable already, has more than one variable group, or has any homogeneous variable groups.
			*/
			TotalDegree(System const& s);


			/**
			Get the random value for start function with index

			\param index The index of the start function for which you want the corresponding random value.
			*/
			template<typename NumT>
			NumT RandomValue(size_t index) const
			{
				return random_values_[index]->Eval<NumT>();
			}


			/**
			Get all the random values, in their Node form.
			*/
			std::vector<std::shared_ptr<node::Rational> > const& RandomValues()
			{
				return random_values_;
			}


			/**
			Get the number of start points for this total degree start system.  This is the Bezout bound for the target system.  Provided here for your convenience.
			*/
			unsigned long long NumStartPoints() const override;

			TotalDegree& operator*=(Nd const& n);

			TotalDegree& operator+=(System const& sys) = delete;
			
			void SanityChecks(System const& s);

		private:

			/**
			Copy the degrees from another system into this one
			*/
			void CopyDegrees(System const& s);

			/**
			Populate the random values of this system.
			*/
			void SeedRandomValues(int num_functions);

			/**
			Generate the functions for this total degree start system.  Assumes the random values, degrees, and variables are already g2g.
			*/
			void GenerateFunctions();


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



