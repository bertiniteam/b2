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
		
		
		
		
		mpfr_float abs()
		{
			return sqrt(abs2());
		}
		
		mpfr_float abs2()
		{
			return pow(real(),2)+pow(imag(),2);
		}
		
		
		void precision(unsigned int prec)
		{
			real_.precision(prec);
			imag_.precision(prec);
		}
		
		unsigned int precision() const
		{
			return real_.precision();
		}
		
	};
	
	
	inline complex operator+(complex lhs, const complex & rhs){
		lhs += rhs;
		return lhs;
	}
	
}


#endif


