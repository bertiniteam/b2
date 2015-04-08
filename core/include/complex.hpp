#ifndef COMPLEX_CLASS_HPP
#define COMPLEX_CLASS_HPP

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>



namespace bertini {
	using boost::multiprecision::mpfr_float;
	
	class mpfr_complex {
		
		
		mpfr_float real_, imag_;
		
	public:
		
		mpfr_complex(mpfr_float re, mpfr_float im) : real_(re), imag_(im)
		{}
		
		mpfr_complex():real_(), imag_(){}
		
		
		inline mpfr_float real() const {return real_;}
		inline mpfr_float imag() const {return imag_;}
		
		void real(mpfr_float new_real){real_ = new_real;}
		void imag(mpfr_float new_imag){imag_ = new_imag;}
		
		
		
		mpfr_float abs()
		{
			return sqrt(abs2());
		}
		
		mpfr_float abs2()
		{
			return pow(real(),2)+pow(imag(),2);
		}
	};
	
	
}


#endif


