//This file is part of Bertini 2.0.
//
//start_system.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//start_system.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with start_system.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  start_system.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015

/**
\file start_system.hpp 

\brief Defines start system types.
*/


#ifndef BERTINI_START_SYSTEM_HPP
#define BERTINI_START_SYSTEM_HPP

#include "system.hpp"
#include "limbo.hpp"


namespace bertini 
{
	namespace start_system{

		/**
		\brief Abstract base class for other start systems.

		Abstract base class for other start systems.  Start systems are special types of systems, to which we know solutions.  We also know how to construct various types of start systems from arbitrary polynomial systems.

		This class provides the empty virtual declarations for necessary override functions for specific start systems, including NumStartPoints (provides an upper bound on the number of solutions to the target system), and the private functions GenerateStartPoint(index), in double and multiple precision.  These two Generate functions are called by the templated non-overridden function StartPoint(index), which calls the appropriate one based on template type.
		*/
		class StartSystem : public System
		{

		public:
			

			virtual mpz_int NumStartPoints() const = 0;

			template<typename T>
			Vec<T> StartPoint(mpz_int index) const
			{
				return GenerateStartPoint(T(),index);
			}

			virtual ~StartSystem() = default;
			
		private:
			virtual Vec<dbl> GenerateStartPoint(dbl,mpz_int index) const = 0;
			virtual Vec<mpfr> GenerateStartPoint(mpfr,mpz_int index) const = 0;

			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<System>(*this);
			}

		};






		/**
		\brief StartSystem for 1-homogeneous polynomial systems.

		The most basic and easy-to-construct start system in Numerical Algebraic Geometry.  

		The total degree start system uses functions of the form \f$x_i^{d_i} - r_i\f$, where \f$i\f$ is the index of the function relative to the system, \f$d_i\f$ is the degree of that function, and \f$r_i\f$ is a random complex number.  This is very similar to the roots of unity, except that they are moved away from being centered around the origin, to being centered around a random complex number.  

		Note that the corresponding target system MUST be square -- have the same number of functions and variables.  The start system cannot be constructed otherwise, particularly because it is written to throw at the moment if not square.

		The start points are accesses by index (mpz_int), instead of being generated all at once.
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
			mpfr RandomValue(size_t index)
			{
				return random_values_[index]->Eval<mpfr>();
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
			mpz_int NumStartPoints() const override;


		private:

			/**
			Get the ith start point, in double precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<dbl> GenerateStartPoint(dbl,mpz_int index) const override;

			/**
			Get the ith start point, in current default precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<mpfr> GenerateStartPoint(mpfr,mpz_int index) const override;

			std::vector<std::shared_ptr<node::Rational> > random_values_; ///< stores the random values for the start functions.  x^d-r, where r is stored in this vector.
			std::vector<mpz_int> degrees_; ///< stores the degrees of the functions.


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


#endif


