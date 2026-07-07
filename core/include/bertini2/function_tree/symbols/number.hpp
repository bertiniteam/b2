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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
//
// silviana amethyst
// University of Wisconsin - Eau Claire
//
//  Created by Collins, James B. on 4/30/15.
//
//
// number.hpp:  Declares the class Number.

/**
\file number.hpp

\brief Provides the Number Node types, including Rational, Complex, and Integer

*/


#ifndef BERTINI_NODE_NUMBER_HPP
#define BERTINI_NODE_NUMBER_HPP


#include "bertini2/function_tree/symbols/symbol.hpp"



namespace bertini {
namespace node{


	/**
	\brief Abstract Number type from which other Numbers derive.

	This class represents constant leaves to a function tree.
	*/
	class Number : public Symbol
	{
	public:

		virtual ~Number() = default;





		

	   


		/**
		\brief Get the degree of this node.

		The degree of a number is always 0.  It's a number.
		*/
		inline
		int Degree(std::shared_ptr<Variable> const& /*v*/ = nullptr) const override
		{
			return 0;
		}

		/**
		\brief Get the degree of this node.

		The degree of a number is always 0.  It's a number.
		*/
		inline
		int Degree(VariableGroup const& /*vars*/) const override
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
		\brief Is this node homogeneous?

		Numbers are always homogeneous
		*/
		bool IsHomogeneous(std::shared_ptr<Variable> const& /*v*/ = nullptr) const override
		{
			return true;
		}
		
		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& /*vars*/) const override
		{
			return true;
		}

		
		/**
		\brief Differentiate a number.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
	protected:

	private:
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<Symbol>(*this);
		}


	};


	/**
	\brief Signed real Integer storage in an expression tree.

	Signed real Integer storage in an expression tree. Consider using a Rational type.
	*/
	class Integer : public Number
	{
	public:
		BERTINI_DEFAULT_VISITABLE()



		/// \brief Defaulted copy constructor.
		Integer(Integer const&) = default;

		~Integer() = default;




		void print(std::ostream & target) const override;

		// negative literals print with a leading '-', so parenthesize like a Negate
		unsigned Precedence() const override
		{
			return true_value_ < 0 ? PrecNegate : PrecAtom;
		}

		/**
		\brief Get the literal value this node represents.
		*/
		mpz_int const& GetValue() const
		{
			return true_value_;
		}

		bool IsLiteralZero() const override
		{
			return true_value_ == 0;
		}

		bool IsLiteralOne() const override
		{
			return true_value_ == 1;
		}

		std::size_t HashImpl() const override
		{
			std::size_t h = typeid(Integer).hash_code();
			HashCombine(h, std::hash<std::string>{}(true_value_.str()));
			return h;
		}
		bool IsSame(Node const& other) const override
		{
			auto o = dynamic_cast<Integer const*>(&other);
			return o && true_value_ == o->true_value_;
		}

		/// \brief Construct (and intern) a Integer node.
		template<typename... Ts>
		static
		std::shared_ptr<Integer> Make(Ts&& ...ts){
			return std::static_pointer_cast<Integer>(Intern(std::shared_ptr<Node>( new Integer(ts...) )));
		}

	private:
		// https://stackoverflow.com/questions/33933550/exception-bad-weak-ptr-while-shared-from-this
		// struct MakeConstructorPublic;
		
		explicit
		Integer(int val) : true_value_(val)
		{}

		explicit
		Integer(mpz_int val) : true_value_(val)
		{}

		explicit
		Integer(std::string const& val) : true_value_(val)
		{}





		mpz_int true_value_;

		friend class boost::serialization::access;

		Integer() = default;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<Number>(*this);
			ar & true_value_;
		}
	};

	// struct Integer::MakeConstructorPublic : public Integer {
	// 	template <typename... Args> 
	// 	MakeConstructorPublic(Args&& ...args) : Integer(std::forward<Args>(args)...)
	// 	{} 
	// };



	/**
	\brief A complex-number literal node in an expression tree.

	Stores an arbitrary-precision **complex** value (complex_mp) -- a real-valued literal is just
	the special case with zero imaginary part.  The value passed at construction time is held as the
	'true value' at its authored precision, and evaluation down- or up-samples from it.  Despite the
	historical "float" name this node carried, it is NOT real-only: it holds a full complex number.
	Prefer a Rational or Integer when the coefficient is exact -- they evaluate faster and to
	arbitrary precision without a stored sample.
	*/
	class Complex : public Number
	{
	public:
		BERTINI_DEFAULT_VISITABLE()



		~Complex() = default;
		




		void print(std::ostream & target) const override;

		// real-valued floats print bare (no complex pair); negative ones get
		// a leading '-', so parenthesize like a Negate.  pairs self-delimit.
		unsigned Precedence() const override
		{
			if (highest_precision_value_.imag() == 0 && highest_precision_value_.real() < 0)
				return PrecNegate;
			return PrecAtom;
		}

		/**
		\brief Get the literal value this node represents, at its stored (highest) precision.
		*/
		complex_mp const& GetValue() const
		{
			return highest_precision_value_;
		}

		bool IsLiteralZero() const override
		{
			return highest_precision_value_.real() == 0 && highest_precision_value_.imag() == 0;
		}

		bool IsLiteralOne() const override
		{
			return highest_precision_value_.real() == 1 && highest_precision_value_.imag() == 0;
		}

		std::size_t HashImpl() const override
		{
			std::size_t h = typeid(Complex).hash_code();
			HashCombine(h, std::hash<std::string>{}(highest_precision_value_.real().str()));
			HashCombine(h, std::hash<std::string>{}(highest_precision_value_.imag().str()));
			return h;
		}
		bool IsSame(Node const& other) const override
		{
			auto o = dynamic_cast<Complex const*>(&other);
			return o
				&& highest_precision_value_.real() == o->highest_precision_value_.real()
				&& highest_precision_value_.imag() == o->highest_precision_value_.imag();
		}

		/// \brief Construct (and intern) a Complex node.
		template<typename... Ts>
		static
		std::shared_ptr<Complex> Make(Ts&& ...ts){
			return std::static_pointer_cast<Complex>(Intern(std::shared_ptr<Node>( new Complex(ts...) )));
		}

	private:

		explicit
		Complex(complex_mp const& val) : highest_precision_value_(val)
		{}

		explicit
		Complex(real_mp const& rval, real_mp const& ival = 0) : highest_precision_value_(rval,ival)
		{}

		explicit
		Complex(std::string const& val) : highest_precision_value_(val)
		{}

		explicit
		Complex(std::string const& rval, std::string const& ival) : highest_precision_value_(rval,ival)
		{}



		complex_mp highest_precision_value_;

		friend class boost::serialization::access;
		Complex() = default;
		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<Number>(*this);
			ar & const_cast<complex_mp &>(highest_precision_value_);
		}
	};


	



	/**
	\brief The Rational number type for Bertini2 expression trees.

	The Rational number type for Bertini2 expression trees.  The "true value" is stored using two mpq_rational numbers from the Boost.Multiprecision library, and the ratio is converted into a double or a complex_mp at evaluate time.
	*/
	class Rational : public Number
	{
	public:
		BERTINI_DEFAULT_VISITABLE()
		
		using mpq_rational = bertini::mpq_rational;  ///< The exact rational type this node stores.

		


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

		// real-valued rationals print bare (no complex pair).  the bare form is
		// textually an expression: a leading '-' parenthesizes like a Negate, and
		// 'p/q' contains a division, so it binds like a Mult (x/(1/3), not x/1/3).
		// complex pairs self-delimit.
		unsigned Precedence() const override
		{
			if (true_value_imag_ == 0)
			{
				if (true_value_real_ < 0)
					return PrecNegate;
				if (denominator(true_value_real_) != 1)
					return PrecMult;
			}
			return PrecAtom;
		}

		/**
		\brief Get the real part of the literal value this node represents.
		*/
		mpq_rational const& GetValueReal() const
		{
			return true_value_real_;
		}

		/**
		\brief Get the imaginary part of the literal value this node represents.
		*/
		mpq_rational const& GetValueImag() const
		{
			return true_value_imag_;
		}

		/**
		\brief Get this exact constant as a number of type NumT, independent of the
		evaluation engine.

		Unlike Eval, this does no caching and never touches the node's stored working
		value --- it is a pure read of the literal.  It matches the literal's conversion:
		double truncation for complex_dbl, and a value at the current thread precision for mpfr.
		*/
		template<typename NumT>
		NumT Value() const
		{
			if constexpr (std::is_same<NumT, complex_dbl>::value)
				return complex_dbl(double(true_value_real_), double(true_value_imag_));
			else
				return NumT(boost::multiprecision::mpfr_float(true_value_real_, ThreadPrecision()),
				            boost::multiprecision::mpfr_float(true_value_imag_, ThreadPrecision()));
		}

		bool IsLiteralZero() const override
		{
			return true_value_real_ == 0 && true_value_imag_ == 0;
		}

		bool IsLiteralOne() const override
		{
			return true_value_real_ == 1 && true_value_imag_ == 0;
		}

		std::size_t HashImpl() const override
		{
			std::size_t h = typeid(Rational).hash_code();
			HashCombine(h, std::hash<std::string>{}(true_value_real_.str()));
			HashCombine(h, std::hash<std::string>{}(true_value_imag_.str()));
			return h;
		}
		bool IsSame(Node const& other) const override
		{
			auto o = dynamic_cast<Rational const*>(&other);
			return o && true_value_real_ == o->true_value_real_ && true_value_imag_ == o->true_value_imag_;
		}



		/// \brief Construct (and intern) a Rational node.
		template<typename... Ts>
		static
		std::shared_ptr<Rational> Make(Ts&& ...ts){
			return std::static_pointer_cast<Rational>(Intern(std::shared_ptr<Node>( new Rational(ts...) )));
		}

	private:

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




		mpq_rational true_value_real_, true_value_imag_;
		Rational() = default;
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<Number>(*this);
			ar & const_cast<mpq_rational &>(true_value_real_);
			ar & const_cast<mpq_rational &>(true_value_imag_);
		}
	};



} // re: namespace node
} // re: namespace bertini

#endif
