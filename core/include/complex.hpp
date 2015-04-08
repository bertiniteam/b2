#ifndef COMPLEX_CLASS_HPP
#define COMPLEX_CLASS_HPP

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>



namespace bertini {
	using boost::multiprecision::mpfr_float;
	
	class complex {
		
		
		mpfr_float real_, imag_;
		
	public:
		
		/**
		 default constructor
		 */
		complex():real_(), imag_(){}
		
		/**
		 two-parameter constructor for building a complex from two high precision numbers
		 */
		complex(mpfr_float re, mpfr_float im) : real_(re), imag_(im)
		{}
		
		
		
		/**
		 two-parameter constructor for building a complex from two strings.
		 */
		complex(std::string re, std::string im) : real_(re), imag_(im)
		{}
		
		
		inline mpfr_float real() const {return real_;}
		inline mpfr_float imag() const {return imag_;}
		
		void real(mpfr_float new_real){real_ = new_real;}
		void imag(mpfr_float new_imag){imag_ = new_imag;}
		
		
		
		complex& operator+=(const complex & rhs)
		{
			real_+=rhs.real_;
			imag_+=rhs.imag_;
			return *this;
		}
		
		
		complex& operator*=(const complex & rhs)
		{
			mpfr_float a, b;
			a = real_*rhs.real_ - imag_*rhs.imag_;
			b = real_*rhs.imag_ + imag_*rhs.real_;
			real_ = a;
			imag_ = b;
			return *this;
		}
		
		
		/**
		 compute the absolute value of the number
		 */
		mpfr_float abs() const
		{
			return sqrt(abs2());
		}
		
		
		/**
		 compute the square of the absolute value of the number
		 */
		mpfr_float abs2() const
		{
			return pow(real(),2)+pow(imag(),2);
		}
		
		
		/**
		 change the precision of this high-precision complex number
		 */
		void precision(unsigned int prec)
		{
			real_.precision(prec);
			imag_.precision(prec);
		}
		
		
		/**
		 get the precision of the high-precision complex number.
		 */
		unsigned int precision() const
		{
			return real_.precision();
		}
		
	};
	
	
	inline mpfr_float abs2(const complex & z)
	{
		return z.abs2();
	}
	
	inline mpfr_float abs(const complex & z)
	{
		return sqrt(abs(z));
	}
	
	
	inline complex operator+(complex lhs, const complex & rhs){
		lhs += rhs;
		return lhs;
	}
	
	inline complex operator*(complex lhs, const complex & rhs){
		lhs *= rhs;
		return lhs;
	}
	
}


#endif


