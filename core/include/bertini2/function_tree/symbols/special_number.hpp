//This file is part of Bertini 2.
//
//special_number.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//special_number.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with special_number.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
// Dani Brake
// University of Notre Dame
//
//  Created by Collins, James B. on 4/30/15.
//
//
// special_number.hpp:  Declares the class SpecialNumber.

/**
\file special_number.hpp

\brief Provides the special numbers \f$\pi\f$, \f$e\f$, 1, 2, 0, as Nodes

*/

#ifndef BERTINI_FUNCTION_TREE_SPECIAL_NUMBER_HPP
#define BERTINI_FUNCTION_TREE_SPECIAL_NUMBER_HPP


#include "bertini2/function_tree/symbols/number.hpp"

#include <cmath>



namespace bertini {
namespace node{
	using ::acos;
	using ::exp;

	
	namespace special_number{

		/**
		\brief The number \f$\pi\f$.

		The number \f$\pi\f$.  Gets its own class because it is such an important number.
		*/
		class Pi : public virtual Number, public virtual NamedSymbol
		{
		public:
			Pi() : NamedSymbol("pi")
			{}

			virtual ~Pi() = default;

		private:
			// Return value of constant
			dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
			
			void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


			mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
			
			void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<NamedSymbol>(*this);
			}
		};



		/**
		\brief The number \f$e\f$.

		The number \f$e\f$.  Gets its own class because it is such an important number.
		*/
		class E : public virtual Number, public virtual NamedSymbol
		{
		public:
			E() : NamedSymbol("e")
			{}

			virtual ~E() = default;



		private:
			// Return value of constant
			dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
			
			void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


			mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
			
			void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<NamedSymbol>(*this);
			}

		};

	} // re: special_number namespace


	/**
	Construct a shared pointer to \f$\pi\f$.
	*/
	std::shared_ptr<Node> Pi();

	/**
	Construct a shared pointer to \f$e\f$.
	*/
	std::shared_ptr<Node> E();

	/**
	Construct a shared pointer to \f$i\f$.
	*/
	std::shared_ptr<Node> I();

	/**
	Construct a shared pointer to \f$2\f$.
	*/
	std::shared_ptr<Node> Two();

	/**
	Construct a shared pointer to \f$1\f$.
	*/
	std::shared_ptr<Node> One();

	/**
	Construct a shared pointer to \f$0\f$.
	*/
	std::shared_ptr<Node> Zero();


} // re: namespace node
} // re: namespace bertini
	
#endif





