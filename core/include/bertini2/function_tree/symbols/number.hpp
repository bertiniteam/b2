//This file is part of Bertini 2.0.
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


#include "function_tree/symbols/symbol.hpp"



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



		void Reset() override
		{
			// nothing to reset here
		}


		

	   


		/**
		Compute the degree with respect to a single variable.

		For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		*/
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return 0;
		}


		int Degree(VariableGroup const& vars) const override
		{
			return 0;
		}

		std::vector<int> MultiDegree(VariableGroup const& vars) const override
		{
			return std::vector<int>(vars.size(), 0);
		}


		void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
		{
			
		}

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
		virtual void precision(unsigned int prec) override
		{
			auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
			val_pair.first.precision(prec);
		}



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
		
		Integer()
		{}

		Integer(int val) : true_value_(val)
		{}

		Integer(mpz_int val) : true_value_(val)
		{}

		explicit
		Integer(std::string const& val) : true_value_(val)
		{}


		~Integer() = default;
		



		void print(std::ostream & target) const override
		{
			target << true_value_;
		}


		/**
		 Differentiates a number.  Should this return the special number Zero?
		 */
		std::shared_ptr<Node> Differentiate() override
		{
			return std::make_shared<Integer>(0);
		}



	private:

		// Return value of constant
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
		{
			return dbl(double(true_value_),0);
		}

		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
		{
			return mpfr(true_value_,0);
		}

		mpz_int true_value_;

		friend class boost::serialization::access;

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
		Float()
		{}


		Float(mpfr const& val) : highest_precision_value_(val)
		{}

		Float(mpfr_float const& rval, mpfr_float const& ival = 0) : highest_precision_value_(rval,ival)
		{}

		Float(std::string const& val) : highest_precision_value_(val)
		{}

		Float(std::string const& rval, std::string const& ival) : highest_precision_value_(rval,ival)
		{}

		~Float() = default;
		




		void print(std::ostream & target) const override
		{
			target << highest_precision_value_;
		}


		/**
		 Differentiates a number.  Should this return the special number Zero?
		 */
		std::shared_ptr<Node> Differentiate() override
		{
			return std::make_shared<Integer>(0);
		}


	private:
		// Return value of constant
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
		{
			return dbl(highest_precision_value_);
		}

		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
		{
			return mpfr(highest_precision_value_);
		}

		mpfr highest_precision_value_;

		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Number>(*this);
			ar & highest_precision_value_;
		}
	};


	



	/**
	\brief The Rational number type for Bertini2 expression trees.

	The Rational number type for Bertini2 expression trees.  The `true value' is stored using two mpq_rational numbers from the Boost.Multiprecision library, and the ratio is converted into a double or a mpfr at evaluate time.
	*/
	class Rational : public virtual Number
	{
	public:


		Rational()
		{}

		Rational(std::string val) : true_value_real_(val), true_value_imag_(0)
		{}

		Rational(std::string val_real, std::string val_imag) : true_value_real_(val_real), true_value_imag_(val_imag)
		{}


		Rational(mpq_rational const& val_real, mpq_rational const& val_imag = 0) : true_value_real_(val_real), true_value_imag_(val_imag)
		{}

		Rational(int, int) = delete;

		~Rational() = default;
		
		static Rational Rand()
		{
			return Rational(RandomRat(),RandomRat());
		}

		static Rational RandReal()
		{
			return Rational(RandomRat());
		}

		void print(std::ostream & target) const override
		{
			target << "(" << true_value_real_ << "," << true_value_imag_ << ")";
		}


		/**
		 Differentiates a number.  
		 */
		std::shared_ptr<Node> Differentiate() override
		{
			return std::make_shared<Integer>(0);
		}



	private:

		// Return value of constant
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
		{
			return dbl(double(true_value_real_),double(true_value_imag_));
		}

		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
		{
			return mpfr(boost::multiprecision::mpfr_float(true_value_real_),boost::multiprecision::mpfr_float(true_value_imag_));
		}

		mpq_rational true_value_real_, true_value_imag_;

		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Number>(*this);
			ar & true_value_real_;
			ar & true_value_imag_;
		}
	};


} // re: namespace node
} // re: namespace bertini

#endif
