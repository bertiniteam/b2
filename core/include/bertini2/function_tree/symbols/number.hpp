//This file is part of Bertini 2.
//
//number.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//number.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with number.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
// number.hpp:  Declares the class Number.

/**
\file number.hpp

\brief Provides the Number Node types, including Rational, Float, and Integer

*/


#ifndef BERTINI_NODE_NUMBER_HPP
#define BERTINI_NODE_NUMBER_HPP


#include "bertini2/function_tree/symbols/symbol.hpp"
#include "bertini2/function_tree/factory.hpp"


namespace bertini {
namespace node{


	/**
	\brief Abstract Number type from which other Numbers derive.

	This class represents constant leaves to a function tree.  FreshEval simply returns
	the value of the constant.
	*/
	class Number : public virtual Symbol
	{
	public:

		virtual ~Number() = default;



		void Reset() const override;


		

	   


		/**
		\brief Get the degree of this node.

		The degree of a number is always 0.  It's a number.
		*/
		inline
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return 0;
		}

		/**
		\brief Get the degree of this node.

		The degree of a number is always 0.  It's a number.
		*/
		inline
		int Degree(VariableGroup const& vars) const override
		{
			return 0;
		}


		/**
		\brief Get the multidegree of this node.

		The degree of a number is always 0.  It's a number.
		*/
		inline
		std::vector<int> MultiDegree(VariableGroup const& vars) const override
		{
			return std::vector<int>(vars.size(), 0);
		}

		/**
		\brief Homogenize this node.

		Homogenization of a number is a trivial operation.  Don't do anything.
		*/
		void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
		{
			
		}

		/**
		\brief Is this node homogeneous?

		Numbers are always homogeneous
		*/
		bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return true;
		}
		
		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& vars) const override
		{
			return true;
		}

		
		/**
		 Change the precision of this variable-precision tree node.
		 
		 \param prec the number of digits to change precision to.
		 */
		void precision(unsigned int prec) const override;

		/**
		\brief Differentiate a number.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
	protected:

	private:
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Symbol>(*this);
		}


	};


	/**
	\brief Signed real Integer storage in an expression tree.

	Signed real Integer storage in an expression tree. Consider using a Rational type.
	*/
	class Integer : public virtual Number
	{
	public:
		
		explicit
		Integer(int val) : true_value_(val)
		{}

		explicit
		Integer(mpz_int val) : true_value_(val)
		{}

		explicit
		Integer(std::string const& val) : true_value_(val)
		{}


		Integer(Integer const&) = default;

		~Integer() = default;
		



		void print(std::ostream & target) const override;


	private:

		// Return value of constant
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


		mpz_int true_value_;

		friend class boost::serialization::access;

		Integer() = default;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Number>(*this);
			ar & true_value_;
		}
	};




	/**
	\brief Number type for storing floating point numbers within an expression tree.  

	 Number type for storing floating point numbers within an expression tree.  The number passed in at construct time is stored as the true value, and evaluation down or up samples from this 'true value'.  Consider using a Rational or Integer if possible.
	*/
	class Float : public virtual Number
	{
	public:

		explicit
		Float(mpfr const& val) : highest_precision_value_(val)
		{}

		explicit
		Float(mpfr_float const& rval, mpfr_float const& ival = 0) : highest_precision_value_(rval,ival)
		{}

		explicit
		Float(std::string const& val) : highest_precision_value_(val)
		{}

		explicit
		Float(std::string const& rval, std::string const& ival) : highest_precision_value_(rval,ival)
		{}

		~Float() = default;
		




		void print(std::ostream & target) const override;





	private:

		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


		mpfr highest_precision_value_;

		friend class boost::serialization::access;
		Float() = default;
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Number>(*this);
			ar & const_cast<mpfr &>(highest_precision_value_);
		}
	};


	



	/**
	\brief The Rational number type for Bertini2 expression trees.

	The Rational number type for Bertini2 expression trees.  The `true value' is stored using two mpq_rational numbers from the Boost.Multiprecision library, and the ratio is converted into a double or a mpfr at evaluate time.
	*/
	class Rational : public virtual Number
	{
	public:

		using mpq_rational = bertini::mpq_rational;

		
		explicit
		Rational(int val) : true_value_real_(val), true_value_imag_(0)
		{}

		explicit
		Rational(int val_real_numerator, int val_real_denomenator,
				 int val_imag_numerator, int val_imag_denomenator) 
					:
					 true_value_real_(val_real_numerator,val_real_denomenator), true_value_imag_(val_imag_numerator,val_imag_denomenator)
		{}

		explicit
		Rational(std::string val) : true_value_real_(val), true_value_imag_(0)
		{}

		explicit
		Rational(std::string val_real, std::string val_imag) : true_value_real_(val_real), true_value_imag_(val_imag)
		{}

		explicit
		Rational(mpq_rational const& val_real, mpq_rational const& val_imag = 0) : true_value_real_(val_real), true_value_imag_(val_imag)
		{}

		Rational(int, int) = delete;

		~Rational() = default;
		
		/**
		\brief Get a random complex rational node.  Numerator and denominator will have about 50 digits.

		\see RandomRat()
		*/
		template<unsigned long Digits = 50>
		static
		Rational Rand()
		{
			return Rational(RandomRat<Digits>(),RandomRat<Digits>());
		}

		/**
		\brief Get a random real rational node.  Numerator and denominator will have about 50 digits.

		\see RandomRat()
		*/
		template<int Digits = 50>
		static
		Rational RandReal()
		{
			return Rational(RandomRat<Digits>(),0);
		}

		void print(std::ostream & target) const override;






	private:

		// Return value of constant
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


		mpq_rational true_value_real_, true_value_imag_;
		Rational() = default;
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Number>(*this);
			ar & const_cast<mpq_rational &>(true_value_real_);
			ar & const_cast<mpq_rational &>(true_value_imag_);
		}
	};


} // re: namespace node
} // re: namespace bertini

#endif
