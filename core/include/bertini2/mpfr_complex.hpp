//This file is part of Bertini 2.
//
//mpfr_complex.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mpfr_complex.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with mpfr_complex.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// danielle brake, university of wisconsin eau claire

/**
\file mpfr_complex.hpp

\brief The main multiprecision complex number type.
*/


#ifndef BERTINI_MPFR_COMPLEX_HPP
#define BERTINI_MPFR_COMPLEX_HPP

#include "bertini2/config.h"



#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include <string>
#include <assert.h>



#include "bertini2/mpfr_extensions.hpp"

namespace bertini{
	using dbl = std::complex<double>;
}


#if USE_BMP_COMPLEX

#include <boost/multiprecision/mpc.hpp>

namespace bertini{

namespace bmp = boost::multiprecision;
using bmp::backends::mpc_complex_backend;

#ifdef BMP_EXPRESSION_TEMPLATES
	using mpfr_complex = bmp::number<mpc_complex_backend<0>, bmp::et_on >;
#else
	using mpfr_complex = bmp::number<mpc_complex_backend<0>, bmp::et_off >;
#endif


	inline auto DefaultPrecision()
	{
		auto p = mpfr_float::default_precision();
		assert(p==mpfr_complex::default_precision() && "precision of real and complex multiprecision numbers have drifted...");
		return p;
	}

	inline void DefaultPrecision(unsigned prec)
	{
		mpfr_float::default_precision(prec);
		mpfr_complex::default_precision(prec);
	}
}






#else // the not-boost-provided complex case, !USE_BMP_COMPLEX

namespace bertini {

	inline auto DefaultPrecision()
	{
		return mpfr_float::default_precision();
	}

	inline void DefaultPrecision(unsigned prec)
	{
		mpfr_float::default_precision(prec);
	}

	/**
	\brief Custom multiple precision complex class.
	
	Custom multiple precision complex class.  Carries arbitrary precision through all defined operations.
	Tested for compatibility with Eigen linear algebra library.

	This class currently uses Boost.Multiprecision -- namely, the mpfr_float type for variable precision.
	This class is serializable using Boost.Serialize.
	
	The precision of a newly-made bertini::custom_complex is whatever current default is, set by DefaultPrecision(...).

	\todo{Implement MPI send/receive commands using Boost.MPI or alternative.}
	*/
	class custom_complex {
		
	private:
		// The real and imaginary parts of the complex number
		mpfr_float real_, imag_;
		
		#ifdef USE_THREAD_LOCAL
			static thread_local mpfr_float temp_[8];
		#else
			static mpfr_float temp_[8];
		#endif

		// Let the boost serialization library have access to the private members of this class.
		friend class boost::serialization::access;
		
		/**
		 \brief Save method for archiving a bertini::custom_complex
		 */
		template<class Archive>
		void save(Archive & ar, const unsigned int version) const
		{
			#ifndef BERTINI_DISABLE_ASSERTS
			assert(real_.precision()==imag_.precision() && "real and imaginary parts at different precision at save time for Boost serialization of bertini::custom_complex");
			#endif

			// note, version is always the latest when saving
			unsigned int temp_precision = real_.precision();
			ar & temp_precision;
			ar & real_;
			ar & imag_;
		}
		
		
		/**
		 \brief Load method for archiving a bertini::custom_complex
		 */
		template<class Archive>
		void load(Archive & ar, const unsigned int version)
		{
			unsigned int temp_precision;
			ar & temp_precision;
			
			this->precision(temp_precision);
			
			ar & real_;
			ar & imag_;
		}
		
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		// this has to be here so that the Boost.Serialization library
		// knows that there are different methods for the serialize() method.
		
		
	public:
		
		/**
		 Default constructor
		 */
		custom_complex() : real_(), imag_(){}
		
		
		
		////////
		//
		//  Construct real-valued complex numbers
		//
		//////////////
		
		/**
		 Single-parameter for constructing a real-valued complex from a single real double number
		 */
		explicit
		custom_complex(double re) : real_(re), imag_(0){}
		
		template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
		custom_complex(T re) : real_(re), imag_(0){}

		template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
		explicit
		custom_complex(T re, T im) : real_(re), imag_(im){}

		custom_complex(mpz_int const& re) : real_(re), imag_(0){}

		explicit
		custom_complex(mpz_int const& re, mpz_int const& im) : real_(re), imag_(im){}

		explicit
		custom_complex(mpq_rational const& re) : real_(re), imag_(0){}

		explicit
		custom_complex(mpq_rational const& re, mpq_rational const& im) : real_(re), imag_(im){}
		/**
		 Single-parameter for constructing a real-valued complex from a single high-precision number
		 */
		
		explicit
		custom_complex(const mpfr_float & re) : real_(re), imag_(0){}


		template<typename T, typename S, typename R, 
		 			typename Q = typename std::enable_if<boost::is_convertible<R, mpfr_float>::value, mpfr_float>::type>
		explicit
		custom_complex(boost::multiprecision::detail::expression<T,S,R> const& other)
		{
			real_ = other;
			imag_ = 0;
		}
		
		/**
		 Single-parameter for constructing a real-valued complex from a convertible single string
		 */
		explicit
		custom_complex(const std::string & re) : real_(re), imag_(0){}
		
		
		
		
		////////
		//
		//  construct complex numbers from two parameters
		//
		//////////////
		
		/**
		 Two-parameter constructor for building a complex from two high precision numbers
		 */
		explicit
		custom_complex(const mpfr_float & re, const mpfr_float & im) : real_(re), imag_(im)
		{}
		
		
		/**
		 Two-parameter constructor for building a complex from two low precision numbers
		 */
		 explicit
		custom_complex(std::complex<double> z) : real_(z.real()), imag_(z.imag())
		{}

		/**
		 Two-parameter constructor for building a complex from two low precision numbers
		 */
		 explicit
		custom_complex(double re, double im) : real_(re), imag_(im)
		{}
		
		
		
		/**
		 Two-parameter constructor for building a complex from two strings.
		 */
		explicit
		custom_complex(const std::string & re, const std::string & im) : real_(re), imag_(im)
		{}
		
		
		/**
		 Mixed two-parameter constructor for building a complex from two strings.
		 */
		 explicit
		custom_complex(const mpfr_float & re, const std::string & im) : real_(re), imag_(im)
		{}
		
		
		/**
		 Mixed two-parameter constructor for building a complex from two strings.
		 */
		 explicit
		custom_complex(const std::string & re, const mpfr_float & im) : real_(re), imag_(im)
		{}
		
		
		
		
		
		
		////////
		//
		//  other constructors
		//
		//////////////
		
		
		/**
		 The move constructor
		 */
		custom_complex(custom_complex && other) : real_(std::move(other.real_)), imag_(std::move(other.imag_))
		{}
		
		
		
		/**
		 The copy constructor
		 */
		custom_complex(const custom_complex & other) : real_(other.real_), imag_(other.imag_)
		{}
		
		/**
		 Enable swapping
		 */
		friend void swap(custom_complex & first, custom_complex& second) noexcept
		{
			using std::swap;
			
			swap(first.real_,second.real_);
			swap(first.imag_,second.imag_);
		}
		
		void swap(custom_complex & other)
		{
			using std::swap;
			swap(*this, other);
		}
		/**
		 Assignment operator
		 */
#if BERTINI_ENABLE_COPY_AND_SWAP
		custom_complex& operator=(custom_complex other)
		{	
			swap(*this,other);
			return *this;
		}
#else
		custom_complex& operator=(custom_complex const& other)
		{	

			if (this != &other)
		    {
				real_ = other.real_;
				imag_ = other.imag_;
			}
			return *this;
		}

		custom_complex& operator=(custom_complex && other) = default;
#endif
		
		template<typename T, typename S, typename R, 
		 			typename Q = typename std::enable_if<boost::is_convertible<R, mpfr_float>::value, mpfr_float>::type>
		custom_complex& operator=(boost::multiprecision::detail::expression<T,S,R> const& other)
		{
			real_ = other;
			imag_ = 0;
			return *this;
		}


		custom_complex& operator=(mpfr_float const& other)
		{
			real_ = other;
			imag_ = 0;
			return *this;
		}

		custom_complex& operator=(mpz_int const& other)
		{
			real_ = other;
			imag_ = 0;
			return *this;
		}
		
		template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
		custom_complex& operator=(T re){real_ = re; imag_ = 0; return *this;}
		
		
		
		////////
		//
		//  getters and setters
		//
		//////////////
		
		/**
		 Get the value of the real part of the complex number
		 */
		inline const mpfr_float& real() const {return real_;}
		
		/**
		 Get the value of the imaginary part of the complex number
		 */
		inline const mpfr_float& imag() const {return imag_;}
		
		/**
		 Set the value of the real part of the complex number
		 */
		inline
		void real(const mpfr_float & new_real){real_ = new_real;}
		
		/**
		 Set the value of the imaginary part of the complex number
		 */
		inline
		void imag(const mpfr_float & new_imag){imag_ = new_imag;}
		
		/**
		 Set the value of the real part of the complex number
		 */
		inline
		void real(int new_real){real_ = new_real;}
		
		/**
		 Set the value of the imaginary part of the complex number
		 */
		inline
		void imag(int new_imag){imag_ = new_imag;}

		/**
		 Set the value of the real part of the complex number
		 */
		inline
		void real(mpz_int new_real){real_ = new_real;}
		
		/**
		 Set the value of the imaginary part of the complex number
		 */
		inline
		void imag(mpz_int new_imag){imag_ = new_imag;}
		
		/**
		 Set the value of the real part of the complex number
		 */
		inline
		void real(mpq_rational new_real){real_ = new_real;}
		
		/**
		 Set the value of the imaginary part of the complex number
		 */
		inline
		void imag(mpq_rational new_imag){imag_ = new_imag;}

		

		/**	
		 Set the value of the real part of the complex number, from a double-quoted string.
		 */
		inline
		void real(const std::string & new_real){real_ = mpfr_float(new_real);}
		
		/**
		 Set the value of the imaginary part of the complex number, from a double-quoted string.
		 */
		inline
		void imag(const std::string & new_imag){imag_ = mpfr_float(new_imag);}
		
		
		
		inline
		void SetZero()
		{
			real_ = 0;
			imag_ = 0;
		}
		
		inline
		void SetOne()
		{
			real_ = 1;
			imag_ = 0;
		}
		
		////////
		//
		//  Some static functions which construct special numbers
		//
		//////////////
		
		
		
		
		/**
		Constuct the number 0.
		*/
		inline static custom_complex zero(){
			return custom_complex(0,0);
		}
		
		/**
		Constuct the number 1.
		*/
		inline static custom_complex one(){
			return custom_complex(1,0);
		}
		
		/**
		Constuct the number 2.
		*/
		inline static custom_complex two(){
			return custom_complex(2,0);
		}
		
		/**
		Constuct the number \f$i\f$.
		*/
		inline static custom_complex i()
		{
			return custom_complex(0,1);
		}
		
		/**
		Constuct the number 0.5.
		*/
		inline static custom_complex half()
		{
			return custom_complex("0.5","0");
		}
		
		/**
		Constuct the number -1.
		*/
		inline static custom_complex minus_one()
		{
			return custom_complex(-1,0);
		}
		
		
		// /**
		//  Produce a random complex number, to default precision.
		//  */
		// inline static custom_complex rand()
		// {
		// 	custom_complex returnme( RandomMp(mpfr_float(-1),mpfr_float(1)), RandomMp(mpfr_float(-1),mpfr_float(1)) );
		// 	returnme /= sqrt( returnme.abs());
		// 	return returnme;
		// }
		
		// inline static custom_complex RandomUnit()
		// {
		// 	custom_complex returnme( RandomMp(mpfr_float(-1),mpfr_float(1)), RandomMp(mpfr_float(-1),mpfr_float(1)) );
		// 	returnme /= returnme.abs();
		// 	return returnme;
		// }
		// /**
		//  Produce a random real number \f$\in [-1,\,1]\f$, to current default precision. 
		//  */
		// inline static custom_complex RandomReal()
		// {
		// 	using std::sqrt;
		// 	custom_complex returnme( RandomMp(mpfr_float(-1),mpfr_float(1)), RandomMp(mpfr_float(-1),mpfr_float(1)) );
		// 	returnme /= sqrt( returnme.abs());
		// 	return returnme;
		// }
		
		
		
		
		
		
		////////
		//
		//  The fundamental arithmetic operators
		//
		//////////////
		
		
		/**
		 Complex addition
		 */
		inline
		custom_complex& operator+=(const custom_complex & rhs)
		{
			real_+=rhs.real_;
			imag_+=rhs.imag_;
			return *this;
		}
		
		inline
		custom_complex& operator+=(const mpz_int & rhs)
		{
			real_+=rhs;
			return *this;
		}

		inline
		custom_complex& operator+=(const mpq_rational & rhs)
		{
			real_+=rhs;
			return *this;
		}

		inline
		custom_complex& operator+=(const mpfr_float & rhs)
		{
			real_+=rhs;
			return *this;
		}


		/**
		 Complex addition, by an integral type.
		 */
		template<typename Int, typename = typename std::enable_if<std::is_integral<Int>::value >::type>
		inline
		custom_complex& operator+=(const Int & rhs)
		{
			real_ += rhs;
			return *this;
		}


		/**
		 Complex subtraction
		 */
		inline
		custom_complex& operator-=(const custom_complex & rhs)
		{
			real_-=rhs.real_;
			imag_-=rhs.imag_;
			return *this;
		}
		
		/**
		 Complex subtraction
		 */
		inline
		custom_complex& operator-=(const mpz_int & rhs)
		{
			real_-=rhs;
			return *this;
		}

		/**
		 Complex subtraction
		 */
		inline
		custom_complex& operator-=(const mpq_rational & rhs)
		{
			real_-=rhs;
			return *this;
		}

		/**
		 Complex subtraction
		 */
		inline
		custom_complex& operator-=(const mpfr_float & rhs)
		{
			real_-=rhs;
			return *this;
		}

		/**
		 Complex subtraction, by an integral type.
		 */

		template<typename Int, typename = typename std::enable_if<std::is_integral<Int>::value >::type>
		inline
		custom_complex& operator-=(const Int & rhs)
		{
			real_ -= rhs;
			return *this;
		}



		/**
		 Complex multiplication.  uses a single temporary variable
		 
		 1 temporary, 4 multiplications
		 */
		inline
		custom_complex& operator*=(const custom_complex & rhs)
		{
			temp_[0].precision(DefaultPrecision());

			temp_[0] = real_*rhs.real_ - imag_*rhs.imag_; // cache the real part of the result
			imag_ = real_*rhs.imag_ + imag_*rhs.real_;
			real_ = temp_[0];
			return *this;
		}
		

		/**
		 Complex multiplication, by an integral type.
		 */
		template<typename Int, typename = typename std::enable_if<std::is_integral<Int>::value >::type>
		inline
		custom_complex& operator*=(const Int & rhs)
		{
			real_ *= rhs;
			imag_ *= rhs;
			return *this;
		}

		template<typename T, typename S, typename R> // ,  typename Q = typename std::enable_if<boost::is_convertible<R, mpfr_float>::value, mpfr_float>::type
		inline
		custom_complex& operator*=(const boost::multiprecision::detail::expression<T,S,R> & rhs)
		{
			real_ *= rhs;
			imag_ *= rhs;
			return *this;
		}
		
		inline
		custom_complex& operator*=(const mpz_int & rhs)
		{
			real_ *= rhs;
			imag_ *= rhs;
			return *this;
		}

		inline
		custom_complex& operator*=(const mpq_rational & rhs)
		{
			real_ *= rhs;
			imag_ *= rhs;
			return *this;
		}

		inline
		custom_complex& operator*=(const mpfr_float & rhs)
		{
			real_ *= rhs;
			imag_ *= rhs;
			return *this;
		}

		/**
		 Complex division.  implemented using two temporary variables
		 */
		inline
		custom_complex& operator/=(const custom_complex & rhs)
		{
			temp_[1].precision(DefaultPrecision());
			temp_[2].precision(DefaultPrecision());

			temp_[1] = rhs.abs2(); // cache the denomenator...
			temp_[2] = real_*rhs.real_ + imag_*rhs.imag_; // cache the numerator of the real part of the result
			imag_ = (imag_*rhs.real_ - real_*rhs.imag_)/temp_[1];
			real_ = temp_[2]/temp_[1];
			
			return *this;
		}
		
		template<typename T, typename S, typename R> // ,  typename Q = typename std::enable_if<boost::is_convertible<R, mpfr_float>::value, mpfr_float>::type
		inline
		custom_complex& operator/=(const boost::multiprecision::detail::expression<T,S,R> & rhs)
		{
			real_ /= rhs;
			imag_ /= rhs;
			
			return *this;
		}

		/**
		 Complex division, by a real mpfr_float.
		 */
		inline
		custom_complex& operator/=(const mpfr_float & rhs)
		{
			real_ /= rhs;
			imag_ /= rhs;
			
			return *this;
		}
		
		/**
		 Complex division, by a real mpz_int.
		 */
		inline
		custom_complex& operator/=(const mpz_int & rhs)
		{
			real_ /= rhs;
			imag_ /= rhs;
			
			return *this;
		}

		/**
		 Complex division, by a real mpq.
		 */
		inline
		custom_complex& operator/=(const mpq_rational & rhs)
		{
			real_ /= rhs;
			imag_ /= rhs;
			
			return *this;
		}


		/**
		 Complex division, by an integral type.
		 */
		template<typename Int, typename = typename std::enable_if<std::is_integral<Int>::value >::type>
		inline
		custom_complex& operator/=(const Int & rhs)
		{
			real_ /= rhs;
			imag_ /= rhs;
			return *this;
		}
		
		/**
		 Complex negation
		 */
		inline
		custom_complex operator-() const
		{
			return custom_complex(-real(), -imag());
		}
		
		
		
		
		
		
		
		
		
		
		
		/**
		 Compute the square of the absolute value of the number
		 */
		inline
		mpfr_float abs2() const
		{
			return real()*real()+imag()*imag();
		}
		
		/**
		 Compute the absolute value of the number
		 */
		inline
		mpfr_float abs() const
		{
			return sqrt(abs2());
		}
		
		
		
		
		
		/**
		 Compute the argument of the complex number, with branch cut according to whatever branch boost chose for their atan2 function.
		 */
		inline
		mpfr_float arg() const
		{
			return boost::multiprecision::atan2(imag(),real());
		}
		
		
		/**
		 Compute the inner product of the number with itself.  this is also the magnitude squared.
		 */
		inline
		mpfr_float norm() const
		{
			return abs2();
		}
		
		/**
		 Compute the complex conjugate of the complex number.
		 */
		inline
		custom_complex conj() const
		{
			return custom_complex(real(), -imag());
		}
		
		
		
		/**
		\brief Is \f$z\f$ a NaN?
		*/
		inline
		bool isnan() const
		{
			using boost::math::isnan;
			if (isnan(real()) || isnan(imag()))
				return true;
			else
				return false;
		}
		
		/**
		\brief Is \f$z\f$ \f$\infty\f$?
		*/
		inline
		bool isinf() const
		{
			using boost::math::isinf;
			using boost::math::isnan;
			if ( (!isnan(real()) && !isnan(imag()))
			    &&
			     ( isinf(real()) ||  isinf(imag()))
			   )
				return true;
			else
				return false;
		}


		/**
		 Change the precision of this high-precision complex number.
		 
		 \param prec the number of digits to change precision to.
		 */
		inline
		void precision(unsigned int prec)
		{
			real_.precision(prec);
			imag_.precision(prec);
		}
		
		
		/**
		 Get the precision of the high-precision complex number.
		 
		 \return the number of digits in the number
		 */
		inline
		unsigned int precision() const
		{
			#ifndef BERTINI_DISABLE_ASSERTS
			assert(real_.precision()==imag_.precision() && "real and imaginary parts at different precision when querying precision.  somehow they got out of sync.");
			#endif

			return real_.precision();
		}
		
		
		
		
		
		
		
		
		
		
		
		
		/**
		 Write a complex number to an output stream.
		 
		 Format complies with the std::complex class -- (re,im).
		 */
		friend std::ostream& operator<<(std::ostream& out, const custom_complex & z)
		{
			out << "(" << z.real() << "," << z.imag() << ")";
			return out;
		}
		
		
		/**
		 Read a complex number from an input stream.
		 
		 Format complies with the std::complex class -- (re,im) or (re) or re.
		 
		 If the read fails because of misplaced parentheses, the stream will be in fail state, and the number will be set to NaN.
		 Function may not tolerate white space in the number.
		 */
		friend std::istream& operator>>(std::istream& in, custom_complex & z)
		{
			std::string gotten;
			in >> gotten;
			
			if (gotten[0]=='(') {
				if (*(gotten.end()-1)!=')') {
					in.setstate(std::ios::failbit);
					z.real("NaN");
					z.imag("NaN");
					return in;
				}
				else{
					// try to find a comma in the string.
					size_t comma_pos = gotten.find(",");
					
					// if the second character, have no numbers in the real part.
					// if the second to last character, have no numbers in the imag part.
					
					if (comma_pos!=std::string::npos){
						if (comma_pos==1 || comma_pos==gotten.size()-2) {
							in.setstate(std::ios::failbit);
							z.real("NaN");
							z.imag("NaN");
							return in;
						}
						else{
							z.real(gotten.substr(1, comma_pos-1));
							z.imag(gotten.substr(comma_pos+1, gotten.size()-2 - (comma_pos)));
							return in;
						}
					}
					// did not find a comma
					else{
						z.real(gotten.substr(1,gotten.size()-2));
						z.imag(0);
						return in;
					}
					
				}
			}
			else{
				z.real(gotten);
				z.imag(0);
				return in;
			}
		}
		
		/**
		Test for exact equality of two complex numbers.  Since they are floating point numbers, this comparison is generally unreliable.
		*/
		inline
		bool operator==(custom_complex const& rhs) const
		{
			return (this->real()==rhs.real()) && (this->imag()==rhs.imag());
		}

		inline
		bool operator!=(custom_complex const& rhs) const
		{
			return !(*this==rhs);
		}

		/**
		When explicitly asked, you can convert a bertini::custom_complex into a std::complex<double>.  But only explicitly.  This conversion is narrowing, and should be avoided.
		*/
		explicit operator std::complex<double> () const
		{
			return std::complex<double>(double(real_), double(imag_));
		}
		
		friend void rand(bertini::custom_complex & a, unsigned num_digits);
		friend void RandomReal(bertini::custom_complex & a, unsigned num_digits);
		friend void RandomComplex(bertini::custom_complex & a, unsigned num_digits);
		friend void RandomUnit(bertini::custom_complex & a, unsigned num_digits);



		friend custom_complex operator/(const mpfr_float & lhs, const custom_complex & rhs);
		// friend custom_complex operator/(const custom_complex & lhs, const mpfr_float & rhs);
		friend custom_complex operator/(const mpz_int & lhs, const custom_complex & rhs);
		// friend custom_complex operator/(const custom_complex & lhs, const mpz_int & rhs);
		friend custom_complex operator/(const mpq_rational & lhs, const custom_complex & rhs);
		// friend custom_complex operator/(const custom_complex & lhs, const mpq_rational & rhs);

		friend custom_complex inverse(const custom_complex & z);

		friend custom_complex exp(const custom_complex & z);
	}; // end declaration of the bertini::custom_complex number class
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 Complex-complex addition.
	 */
	inline custom_complex operator+(custom_complex lhs, const custom_complex & rhs){
		lhs += rhs;
		return lhs;
	}
	
	/**
	 Complex-real addition.
	 */
	inline custom_complex operator+(custom_complex lhs, const mpfr_float & rhs)
	{
		lhs.real(lhs.real()+rhs);
		return lhs;
	}
	
	/**
	 Real-complex addition.
	 */
	inline custom_complex operator+(const mpfr_float & lhs, custom_complex rhs)
	{
		return rhs+lhs;
	}

	/**
	 Complex-real addition.
	 */
	inline custom_complex operator+(custom_complex lhs, const mpz_int & rhs)
	{
		lhs.real(lhs.real()+rhs);
		return lhs;
	}
	
	/**
	 Real-complex addition.
	 */
	inline custom_complex operator+(const mpz_int & lhs, custom_complex rhs)
	{
		return rhs+lhs;
	}
	
	/**
	 Complex-real addition.
	 */
	template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
	inline custom_complex operator+(custom_complex lhs, T const& rhs)
	{
		lhs += rhs;
		return lhs;
	}
	
	/**
	 Real-complex addition.
	 */
	template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
	inline custom_complex operator+(T const& lhs, custom_complex rhs)
	{
		rhs += lhs;
		return rhs;
	}
	
	
	
	
	/**
	 Complex-complex subtraction
	 */
	inline custom_complex operator-(custom_complex lhs, const custom_complex & rhs){
		lhs -= rhs;
		return lhs;
	}
	
	/**
	 Complex-real subtraction
	 */
	inline custom_complex operator-(custom_complex lhs, const mpfr_float & rhs)
	{
		lhs -= rhs;
		return lhs;
	}
	
	/**
	 Real-complex subtraction
	 */
	inline custom_complex operator-(const mpfr_float & lhs, custom_complex rhs)
	{
		rhs -= lhs;
		return -rhs;
	}
	
	/**
	 Complex-real subtraction
	 */
	inline custom_complex operator-(custom_complex lhs, const mpz_int & rhs)
	{
		lhs -= rhs;
		return lhs;
	}
	
	/**
	 Real-complex subtraction
	 */
	inline custom_complex operator-(const mpz_int & lhs, custom_complex rhs)
	{
		rhs -= lhs;
		return -rhs;
	}

	/**
	 Complex-integer subtraction
	 */
	template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
	inline custom_complex operator-(custom_complex lhs, T rhs)
	{
		lhs -= rhs;
		return lhs;
	}
	
	/**
	 Integer-complex subtraction
	 */
	template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
	inline custom_complex operator-(T lhs, custom_complex rhs)
	{
		rhs -= lhs;
		return -rhs;
	}
	


	
	/**
	 Complex-complex multiplication
	 */
	inline custom_complex operator*(custom_complex lhs, const custom_complex & rhs){
		lhs *= rhs;
		return lhs;
	}
	
	/**
	 Complex-real multiplication
	 */
	inline custom_complex operator*(custom_complex lhs, const mpfr_float & rhs)
	{
		lhs.real(lhs.real()*rhs);
		lhs.imag(lhs.imag()*rhs);
		return lhs;
	}
	
	/**
	 Real-complex multiplication
	 */
	inline custom_complex operator*(const mpfr_float & lhs, custom_complex rhs)
	{
		return rhs*lhs; // it commutes!
	}
	
	/**
	 Complex-integer multiplication
	 */
	inline custom_complex operator*(custom_complex lhs, const mpz_int & rhs)
	{
		lhs.real(lhs.real()*rhs);
		lhs.imag(lhs.imag()*rhs);
		return lhs;
	}
	
	/**
	 Integer-complex multiplication
	 */
	inline custom_complex operator*(const mpz_int & lhs, custom_complex rhs)
	{
		return rhs*lhs; // it commutes!
	}
	
	
	/**
	 Complex-integer multiplication
	 */
	template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
	inline custom_complex operator*(custom_complex lhs, T const& rhs)
	{
		lhs *= rhs;
		return lhs;
	}
	
	/**
	 Integer-complex multiplication
	 */
	template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
	inline custom_complex operator*(T const& lhs, custom_complex rhs)
	{
		rhs *= lhs;
		return rhs; // it commutes!
	}

	template<typename T, typename S, typename R, 
		 			typename Q = typename std::enable_if<boost::is_convertible<R, mpfr_float>::value, mpfr_float>::type>
	inline custom_complex operator*(custom_complex lhs, const boost::multiprecision::detail::expression<T,S,R> & rhs)
	{
		lhs*=rhs;
		return lhs;
	}

	template<typename T, typename S, typename R, 
		 			typename Q = typename std::enable_if<boost::is_convertible<R, mpfr_float>::value, mpfr_float>::type>
	inline custom_complex operator*(const boost::multiprecision::detail::expression<T,S,R> & lhs, custom_complex rhs)
	{
		rhs*=lhs;
		return rhs;
	}
	
	/**
	 Complex-complex division
	 */
	inline custom_complex operator/(custom_complex lhs, const custom_complex & rhs){
		lhs /= rhs;
		return lhs;
	}
	
	/**
	 Real-complex division
	 */
	inline custom_complex operator/(const mpfr_float & lhs, const custom_complex & rhs)
	{

		custom_complex::temp_[3].precision(DefaultPrecision());
		custom_complex::temp_[3] = rhs.abs2();
		return custom_complex(lhs*rhs.real()/custom_complex::temp_[3], -lhs*rhs.imag()/custom_complex::temp_[3]);
	}

	/**
	 Complex-real division
	 */
	inline custom_complex operator/(custom_complex lhs, const mpfr_float & rhs)
	{
		lhs /= rhs;
		return lhs;
	}
	
	



	/**
	 Integer-complex division
	 */
	inline custom_complex operator/(const mpz_int & lhs, const custom_complex & rhs)
	{
		custom_complex::temp_[4].precision(DefaultPrecision());
		custom_complex::temp_[4] = rhs.abs2();
		return custom_complex(lhs*rhs.real()/custom_complex::temp_[4], -lhs*rhs.imag()/custom_complex::temp_[4]);
	}
	
	/**
	 Complex-integer division
	 */
	inline custom_complex operator/(custom_complex lhs, const mpz_int & rhs)
	{
		lhs/=rhs;
		return lhs;
	}
	
	/**
	 Rational-complex division
	 */
	inline custom_complex operator/(const mpq_rational & lhs, const custom_complex & rhs)
	{
		custom_complex::temp_[4].precision(DefaultPrecision());
		custom_complex::temp_[4] = rhs.abs2();
		return custom_complex(lhs*rhs.real()/custom_complex::temp_[4], -lhs*rhs.imag()/custom_complex::temp_[4]);
	}
	
	/**
	 Complex-Rational division
	 */
	inline custom_complex operator/(custom_complex lhs, const mpq_rational & rhs)
	{
		lhs/=rhs;
		return lhs;
	}

	/**
	 Integer-complex division
	 */
	template<typename T, typename  = typename std::enable_if<std::is_integral<T>::value >::type>
	inline custom_complex operator/(T const& lhs, const custom_complex & rhs)
	{
		static 
#ifdef USE_THREAD_LOCAL
		thread_local 
#endif
		mpfr_float temp;


		temp.precision(DefaultPrecision()); // this precision change disconcerts me.  does it interrupt the cache pipeline a bunch to do it?  
		
		temp = rhs.abs2();
		return custom_complex(lhs*rhs.real()/temp, -lhs*rhs.imag()/temp);
	}
	
	/**
	 Complex-integer division
	 */
	template<typename T, typename = typename std::enable_if<std::is_integral<T>::value >::type>
	inline custom_complex operator/(custom_complex lhs, T const& rhs)
	{
		lhs/=rhs;
		return lhs;
	}


	








	/**
	 Get the real part of a complex number
	 */
	inline const mpfr_float& real(const custom_complex & z)
	{
		return z.real();
	}
	
	/**
	 Get the imaginary part of a complex number
	 */
	inline const mpfr_float& imag(const custom_complex & z)
	{
		return z.imag();
	}
	
	
	/**
	 Conjugate a complex number
	 */
	inline custom_complex conj(const custom_complex & z)
	{
		return z.conj();
	}
	
	/**
	 \brief The C++ norm of complex number.

	 Mathematically we think of this as the square of the absolute value.
	*/	
	inline mpfr_float norm(const custom_complex & z)
	{
		return z.norm();
	}

	
	/**
	 Compute the square of the absolute value of a complex number
	 */
	inline mpfr_float abs2(const custom_complex & z)
	{
		return z.abs2();
	}
	
	/**
	 Compute the absolute value of a complex number.
	 */
	inline mpfr_float abs(const custom_complex & z)
	{
		return boost::multiprecision::sqrt(abs2(z));
	}
	
	
	/**
	 Compute the argument of a complex number, with branch cut determined by the  atan2  function.
	 */
	inline mpfr_float arg(const custom_complex & z)
	{
		return boost::multiprecision::atan2(z.imag(),z.real());
	}
	
	
	
	
	/**
	 Compute the inverse of a complex number
	 */
	inline custom_complex inverse(const custom_complex & z)
	{
		custom_complex::temp_[6].precision(DefaultPrecision());
		custom_complex::temp_[6] = z.abs2();
		
		return custom_complex(z.real()/custom_complex::temp_[6], -z.imag()/custom_complex::temp_[6]);
	}
	
	
	/**
	 Compute the square of a complex number
	 
	 4 multiplications
	 1 creation of a mpfr_float
	 */
	inline custom_complex square(const custom_complex & z)
	{
		return custom_complex(z.real()*z.real() - z.imag()*z.imag(), mpfr_float(2)*z.real()*z.imag());
	}
	
	
	
	/**
	 Compute the cube of a complex number
	 
	 10 multiplications
	 2 creations of an mpfr_float.
	 
	 This could use fewer multiplications if it used more temporaries
	 */
	inline custom_complex cube(const custom_complex & z)
	{
		//		return custom_complex(x^3 - 3*x*y^2, 3*x^2*y - y^3); // this deliberately left in for the equation.
		return custom_complex(pow(z.real(),3) - 3*z.real()*pow(z.imag(),2),
					   3*pow(z.real(),2)*z.imag() - pow(z.imag(),3));
		
	}
	
	
	
	/**
	 Compute +,- integral powers of a complex number.

	 This function recursively calls itself if the power is negative, by computing the power on the inverse.
	 */
	inline custom_complex pow(const custom_complex & z, int power)
	{
		if (power < 0) {
			return pow(inverse(z), -power);
		}
		else if (power==0)
			return custom_complex(1,0);
		else if(power==1)
			return z;
		else if(power==2)
			return z*z;
		else if(power==3)
			return z*z*z;
		else
		{
			unsigned int p(power);
			custom_complex result(1,0), z_to_the_current_power_of_two = z;
			// have copy of p in memory, can freely modify it.
			do {
				if ( (p & 1) == 1 ) { // get the lowest bit of the number
					result *= z_to_the_current_power_of_two;
				}
				z_to_the_current_power_of_two *= z_to_the_current_power_of_two; // square z_to_the_current_power_of_two
			} while (p  >>= 1);
			
			return result;
		}
	}
	
	
	
	/**
	 Construct a complex number from magnitude and angle.
	 */
	inline custom_complex polar(const mpfr_float & rho, const mpfr_float & theta)
	{
		return custom_complex(rho*cos(theta), rho*sin(theta));
	}
	
	
	/**
	 Compute the square root of a complex number, using branch cut along the -x axis.
	 */
	inline custom_complex sqrt(const custom_complex & z)
	{
		return bertini::polar(sqrt(abs(z)), arg(z)/2);
	}
	
	
	/**
	 Compute e^z for complex z.
	 */
	inline custom_complex exp(const custom_complex & z)
	{
		custom_complex::temp_[7].precision(DefaultPrecision());
		custom_complex::temp_[7] = exp(real(z));
		return custom_complex(custom_complex::temp_[7] * cos(imag(z)), custom_complex::temp_[7] * sin(imag(z)));
	}
	
	/**
	 Compute sine of a complex number
	 */
	inline custom_complex sin(const custom_complex & z)
	{
		return (exp(custom_complex::i()*z) - exp(-custom_complex::i()*z)) / custom_complex::i() / 2;
	}
	
	/**
	 Compute cosine of a complex number
	 */
	inline custom_complex cos(const custom_complex & z)
	{
		return (exp(custom_complex::i()*z) + exp(-custom_complex::i()*z))  / 2;
	}
	
	/**
	 Compute tangent of a complex number
	 */
	inline custom_complex tan(const custom_complex & z)
	{
		return sin(z) / cos(z);
	}
	
	
	/**
	 Compute hyperbolic sine of a complex number
	 */
	inline custom_complex sinh(const custom_complex & z)
	{
		return (exp(z) - exp(-z)) / 2;
	}
	
	/**
	 Compute hyperbolic cosine of a complex number
	 */
	inline custom_complex cosh(const custom_complex & z)
	{
		return (exp(z) + exp(-z))  / 2;
	}
	
	/**
	 Compute hyperbolic tangent of a complex number
	 */
	inline custom_complex tanh(const custom_complex & z)
	{
		return (sinh(z) / cosh(z));
	}
	
	
	
	/**
	 Complex logarithm base e.
	 */
	inline custom_complex log(const custom_complex & z)
	{
		return custom_complex( log(abs(z)), arg(z));
	}
	

	
	/**
	 Compute c^z, for c,z complex numbers
	 */
	inline custom_complex pow(const custom_complex & z, const custom_complex & c)
	{
		return exp(c * log(z));
	}
	
	/**
	 Compute c^z, for c,z complex numbers
	 */
	inline custom_complex pow(const custom_complex & z, const mpfr_float & c)
	{
		return exp(c * log(z));
	}

	template<typename T, typename S, typename R, 
		 			typename Q = typename std::enable_if<boost::is_convertible<R, mpfr_float>::value, mpfr_float>::type>
	inline custom_complex pow(const custom_complex & z, const boost::multiprecision::detail::expression<T,S,R> & c)
	{
		return exp(c * log(z));
	}

	
	/**
	 Inverse sine of complex number
	 */
	inline custom_complex asin(const custom_complex & z)
	{
		return (-custom_complex::i()) * log( custom_complex::i()*z + sqrt( 1 - pow(z,2)) );
	}
	
	
	/**
	 Inverse cosine of complex number
	 */
	inline custom_complex acos(const custom_complex & z)
	{
		return -custom_complex::i() * log( z + custom_complex::i()*sqrt( 1 - pow(z,2) ) );
	}
	
	
	
	
	/**
	 Inverse tangent of complex number
	 */
	inline custom_complex atan(const custom_complex & z)
	{
		return custom_complex::i()/2 * log( (custom_complex::i() + z) / (custom_complex::i() - z) );
	}
	
	
	
	
	/**
	 Inverse hyperbolic sine of complex number
	 */
	inline custom_complex asinh(const custom_complex & z)
	{
		return log( z + sqrt( square(z)+1 )  );
	}
	
	/**
	 Inverse hyperbolic cosine of complex number
	 */
	inline custom_complex acosh(const custom_complex & z)
	{
		return log(  z + sqrt( square(z)-1 ) );
	}
	
	/**
	 Inverse hyperbolic tangent of complex number
	 */
	inline custom_complex atanh(const custom_complex & z)
	{
		return log( (1+z)/(1-z) )/2;
	}
	
	
	
	using mpfr_complex = bertini::custom_complex;

} // namespace bertini








namespace bertini{

	/** 
	\brief Get the precision of a number.

	For mpfr_floats, this calls the precision member method for mpfr_float.
	*/
	inline
	auto Precision(mpfr_complex const& num)
	{
		return num.precision();
	}


	/** 
	\brief Change the precision of a number.

	For mpfr_floats, this calls the precision member method for mpfr_float.
	*/
	inline void Precision(mpfr_complex & num, unsigned prec)
	{
		num.precision(prec);
	}

	inline
	bool isnan(mpfr_complex const& num){return isnan(num.real()) || isnan(num.imag());};




// moved to random.hpp.  delete this commented-out code
	// inline 
	// void RandomReal(bertini::mpfr_complex & a, unsigned num_digits)
	// {
	// 	a.precision(num_digits);
	// 	RandomMp(a.real_,num_digits);
	// 	a.imag_ = 0;
	// }

	// inline 
	// void rand(bertini::mpfr_complex & a, unsigned num_digits)
	// {
	// 	a.precision(num_digits);
	// 	RandomMp(a.real_,num_digits);
	// 	RandomMp(a.imag_,num_digits);
	// }

	// inline 
	// void RandomComplex(bertini::mpfr_complex & a, unsigned num_digits)
	// {
	// 	rand(a,num_digits);
	// }

	// inline 
	// bertini::mpfr_complex RandomComplex(unsigned num_digits)
	// {
	// 	bertini::mpfr_complex z;
	// 	RandomComplex(z, num_digits);
	// 	return z;
	// }

	// inline 
	// void RandomUnit(bertini::mpfr_complex & a, unsigned num_digits)
	// {
	// 	auto prev_precision = DefaultPrecision();

	// 	a.precision(num_digits);
	// 	RandomMp(a.real_,num_digits);
	// 	RandomMp(a.imag_,num_digits);
	// 	a /= abs(a);

	// 	DefaultPrecision(prev_precision);
	// }

	// inline 
	// bertini::mpfr_complex RandomUnit(unsigned num_digits)
	// {
	// 	bertini::mpfr_complex a;
	// 	RandomUnit(a,num_digits);
	// 	return a;
	// }


	
} // re: namespace bertini






#endif // temporary to disable the entire interior of the file




#endif




