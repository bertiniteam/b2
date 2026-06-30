//This file is part of Bertini 2.
//
//total_degree_binomial.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//total_degree_binomial.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with total_degree_binomial.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file total_degree_binomial.hpp

\brief Defines the TotalDegreeBinomial start system type.
*/

#pragma once

#include "bertini2/system/start_base.hpp"
#include "bertini2/system/start/utility.hpp"

namespace bertini
{
	namespace start_system{


		/**
		\brief Roots-of-unity start system for 1-homogeneous polynomial systems.

		A basic and easy-to-construct start system in Numerical Algebraic Geometry.

		This start system uses functions of the form \f$x_i^{d_i} - r_i\f$, where \f$i\f$ is the index of the function relative to the system, \f$d_i\f$ is the degree of that function, and \f$r_i\f$ is a random complex number.  Its start points are \f$r_i^{1/d_i}\f$ times the \f$d_i\f$-th roots of unity -- a structured (roots-of-unity) grid, only the overall scale/phase per variable being randomized.  For start points in *generic* position (random linear products), use TotalDegreeLinearProduct instead.

		Note that the corresponding target system MUST be square -- have the same number of functions and variables.  The start system cannot be constructed otherwise, particularly because it is written to throw at the moment if not square.

		The start points are accesses by index (unsigned long long), instead of being generated all at once.
		*/
		class TotalDegreeBinomial : public StartSystem
		{
		public:
			TotalDegreeBinomial() = default;
			virtual ~TotalDegreeBinomial() = default;

			/**
			 Constructor for making a roots-of-unity start system from a polynomial system

			 \throws std::runtime_error, if the input target system is not square, is not polynomial, has a path variable already, has more than one variable group, or has any homogeneous variable groups.
			*/
			TotalDegreeBinomial(System const& s);


			/**
			Get the random value for start function with index

			\param index The index of the start function for which you want the corresponding random value.
			*/
			template<typename NumT>
			NumT RandomValue(size_t index) const
			{
				// A direct read of the literal constant --- no node-level evaluation.  The Complex
				// node stores a complex_mp at its (max) creation precision; convert to NumT.
				auto const& v = random_values_[index]->GetValue();
				if constexpr (std::is_same<NumT, complex_dbl>::value)
					return complex_dbl(double(v.real()), double(v.imag()));
				else
					return NumT(v);
			}


			/**
			Get all the random values, in their Node form.
			*/
			std::vector<std::shared_ptr<node::Complex> > const& RandomValues()
			{
				return random_values_;
			}


			/**
			Get the number of start points for this roots-of-unity start system.  This is the Bezout bound for the target system.  Provided here for your convenience.
			*/
			unsigned long long NumStartPoints() const override;

			TotalDegreeBinomial& operator*=(Nd const& n);

			TotalDegreeBinomial& operator+=(System const& sys) = delete;

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
			Generate the functions for this roots-of-unity start system.  Assumes the random values, degrees, and variables are already g2g.
			*/
			void GenerateFunctions();


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

			std::vector<std::shared_ptr<node::Complex> > random_values_; ///< stores the random values for the start functions.  x^d-r, where r is stored in this vector (a literal complex_mp, modulus near 1).
			std::vector<unsigned long long> degrees_; ///< stores the degrees of the functions.


			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned /*version*/) {
				ar & boost::serialization::base_object<StartSystem>(*this);
				ar & random_values_;
				ar & degrees_;
			}

		};
	}
}
