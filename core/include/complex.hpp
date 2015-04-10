#ifndef COMPLEX_CLASS_HPP
#define COMPLEX_CLASS_HPP

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>

#include <assert.h>

namespace bertini {
	using boost::multiprecision::mpfr_float;
	
	class complex {
		
		// the real and imaginary parts of the complex number
		mpfr_float real_, imag_;
		
	public:
		
		/**
		 default constructor
		 */
		complex():real_(), imag_(){}
		
		/**
		 two-parameter constructor for building a complex from two high precision numbers
		 */
		complex(const mpfr_float & re, const mpfr_float & im) : real_(re), imag_(im)
		{}
		
		
		
		/**
		 two-parameter constructor for building a complex from two strings.
		 */
		complex(std::string re, const std::string & im) : real_(re), imag_(im)
		{}
		
		
		/**
		 mixed two-parameter constructor for building a complex from two strings.
		 */
		complex(const mpfr_float & re, const std::string & im) : real_(re), imag_(im)
		{}
		
		
		/**
		 mixed two-parameter constructor for building a complex from two strings.
		 */
		complex(const std::string & re, const mpfr_float & im) : real_(re), imag_(im)
		{}
		
		
		
		/**
		 the move constructor
		 */
		complex(complex&& other) : complex() // initialize via default constructor, C++11 only
		{
			swap(*this, other);
		}
		
		
		
		/**
		 the copy constructor
		 */
		complex(const complex & other) : real_(other.real_), imag_(other.imag_)
		{}
		
		/**
		 enable swapping
		 */
		friend void swap(complex& first, complex& second) // nothrow
		{
			using std::swap;
			
			swap(first.real_,second.real_);
			swap(first.imag_,second.imag_);
		}
		
		/**
		 assignment operator
		 */
		complex& operator=(complex other)
		{
			swap(*this, other);
			return *this;
		}
		
		/**
		 get the value of the real part of the complex number
		 */
		inline mpfr_float real() const {return real_;}
		
		/**
		 get the value of the imaginary part of the complex number
		 */
		inline mpfr_float imag() const {return imag_;}
		
		/**
		 Set the value of the real part of the complex number
		 */
		void real(const mpfr_float & new_real){real_ = new_real;}
		
		/**
		 Set the value of the imaginary part of the complex number
		 */
		void imag(const mpfr_float & new_imag){imag_ = new_imag;}
		
		/**
		 Set the value of the real part of the complex number, from a double-quoted string.
		 */
		void real(const std::string & new_real){real_ = mpfr_float(new_real);}
		
		/**
		 Set the value of the imaginary part of the complex number, from a double-quoted string.
		 */
		void imag(const std::string & new_imag){imag_ = mpfr_float(new_imag);}
		
		
		
		
		
		
		
		inline static complex i()
		{
			return complex("0.0","1.0");
		}
		
		inline static complex half()
		{
			return complex("0.5","0.0");
		}
		
		inline static complex minus_one()
		{
			return complex("-1.0","0.0");
		}
		
		inline static complex zero(){
			return complex("0.0","0.0");
		}
		
		inline static complex one(){
			return complex("1.0","0.0");
		}
		
		inline static complex two(){
			return complex("2.0","0.0");
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		/**
		 complex addition
		 */
		complex& operator+=(const complex & rhs)
		{
			real_+=rhs.real_;
			imag_+=rhs.imag_;
			return *this;
		}
		
		
		/**
		 complex subtraction
		 */
		complex& operator-=(const complex & rhs)
		{
			real_-=rhs.real_;
			imag_-=rhs.imag_;
			return *this;
		}
		
		
		/**
		 complex multiplication.  uses a single temporary variable
		 
		 1 temporary, 4 multiplications
		 */
		complex& operator*=(const complex & rhs)
		{
			mpfr_float a = real_*rhs.real_ - imag_*rhs.imag_; // cache the real part of the result
			imag_ = real_*rhs.imag_ + imag_*rhs.real_;
			real_ = a;
			return *this;
		}
		
		
		/**
		 complex division.  implemented using two temporary variables
		 */
		complex& operator/=(const complex & rhs)
		{
			mpfr_float d = rhs.abs2();
			mpfr_float a = real_*rhs.real_ + imag_*rhs.imag_; // cache the numerator of the real part of the result
			imag_ = imag_*rhs.real_ - real_*rhs.imag_/d;
			real_ = a/d;
			
			return *this;
		}
		
		
		
		
		
		
		
		
		/**
		 write a complex number to an output stream.
		 
		 format complies with the std::complex class -- (re,im).
		 */
		friend std::ostream& operator<<(std::ostream& out, const complex & z)
		{
			out << "(" << z.real() << "," << z.imag() << ")";
			return out;
		}
		
		
		/**
		 read a complex number from an input stream.
		 
		 format complies with the std::complex class -- (re,im) or (re) or re.
		 
		 if the read fails because of misplaced parentheses, the stream will be in fail state, and the number will be set to NaN.
		 function does NOT tolerate white space in the number.
		 */
		friend std::istream& operator>>(std::istream& in, complex & z)
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
						z.imag("0.0");
						return in;
					}
					
				}
			}
			else{
				z.real(gotten);
				z.imag("0.0");
				return in;
			}
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
		 get the argument of the complex number, with branch cut according to whatever branch boost chose for their atan2 function.
		 */
		mpfr_float arg() const
		{
			return boost::multiprecision::atan2(imag(),real());
		}
		
		
		/**
		 get the inner product of the number with itself.
		 */
		mpfr_float norm() const
		{
			return abs2();
		}
		
		/**
		 get the complex conjugate of the complex number.
		 */
		complex conj() const
		{
			return complex(real(), -imag());
		}
		
		
		
		
		
		
		/**
		 change the precision of this high-precision complex number. 
		 
		 \param prec the number of digits to change precision to.
		 */
		void precision(unsigned int prec)
		{
			real_.precision(prec);
			imag_.precision(prec);
		}
		
		
		/**
		 get the precision of the high-precision complex number.
		 
		 \return the number of digits in the number
		 */
		unsigned int precision() const
		{
			assert(real_.precision()==imag_.precision());
			
			return real_.precision();
		}
		
	};
	
	
	/**
	 compute the square of the absolute value of a complex number
	 */
	inline mpfr_float abs2(const complex & z)
	{
		return z.abs2();
	}
	
	/**
	 compute the absolute value of a complex number.
	 */
	inline mpfr_float abs(const complex & z)
	{
		return boost::multiprecision::sqrt(abs2(z));
	}
	
	
	/**
	 compute the argument of a complex number, with branch cut determined by the  atan2  function.
	 */
	inline mpfr_float arg(const complex & z)
	{
		return boost::multiprecision::atan2(z.imag(),z.real());
	}
	
	
	
	
	
	
	
	
	/**
	 complex-complex multiplication
	 */
	inline complex operator+(complex lhs, const complex & rhs){
		lhs += rhs;
		return lhs;
	}
	
	/**
	 complex-real multiplication
	 */
	inline complex operator+(complex lhs, const mpfr_float & rhs)
	{
		lhs.real(lhs.real()+rhs);
		return lhs;
	}
	
	/**
	 real-complex multiplication
	 */
	inline complex operator+(const mpfr_float & lhs, complex rhs)
	{
		return rhs+lhs;
	}
	
	
	
	
	
	/**
	 complex-complex subtraction
	 */
	inline complex operator-(complex lhs, const complex & rhs){
		lhs -= rhs;
		return lhs;
	}
	
	/**
	 complex-real subtraction
	 */
	inline complex operator-(complex lhs, const mpfr_float & rhs)
	{
		lhs.real(lhs.real()-rhs);
		return lhs;
	}
	
	/**
	 real-complex subtraction
	 */
	inline complex operator-(const mpfr_float & lhs, complex rhs)
	{
		return rhs-lhs;
	}
	
	
	
	
	/**
	   complex-complex multiplication
	   */
	inline complex operator*(complex lhs, const complex & rhs){
		lhs *= rhs;
		return lhs;
	}
	
	/**
	 complex-real multiplication
	 */
	inline complex operator*(complex lhs, const mpfr_float & rhs)
	{
		lhs.real(lhs.real()*rhs);
		lhs.imag(lhs.imag()*rhs);
		return lhs;
	}
	
	/**
	 real-complex multiplication
	 */
	inline complex operator*(const mpfr_float & lhs, complex rhs)
	{
		return rhs*lhs;
	}
	
	
	
	
	
	/**
	 complex-complex division
	 */
	inline complex operator/(complex lhs, const complex & rhs){
		lhs /= rhs;
		return lhs;
	}
	
	/**
	 complex-real division
	 */
	inline complex operator/(complex lhs, const mpfr_float & rhs)
	{
		lhs.real(lhs.real()/rhs);
		lhs.imag(lhs.imag()/rhs);
		return lhs;
	}
	
	/**
	 real-complex division
	 */
	inline complex operator/(const mpfr_float & lhs, const complex & rhs)
	{
		mpfr_float d = rhs.abs2();
		return complex(lhs*rhs.real()/d, -lhs*rhs.imag()/d);
	}
	
	
	
	
	
	
	/**
	 compute the inverse of a complex number
	 */
	inline complex inverse(const complex & z)
	{
		mpfr_float d =z.abs2();
		
		return complex(z.real()/d, -z.imag()/d);
	}
	
	
	/**
	 compute the square of a complex number
	 
	 4 multiplications
	 1 creation of a mpfr_float
	 */
	inline complex square(const complex & z)
	{
		return complex(z.real()*z.real() - z.imag()*z.imag(), mpfr_float("2.0")*z.real()*z.imag());
	}
	
	/**
	 compute the cube of a complex number
	 
	 10 multiplications
	 2 creations of an mpfr_float.
	 */
	inline complex cube(const complex & z)
	{
		//		return complex(x^3 - 3*x*y^2, 3*x^2*y - y^3); // this deliberately left in for the equation.
		return complex(pow(z.real(),3) - mpfr_float("3.0")*z.real()*pow(z.imag(),2),
					   mpfr_float("3.0")*pow(z.real(),2)*z.imag() - pow(z.imag(),3));
		
	}
	
	
	
	/**
	 compute +- integral powers of a complex number.
	 */
	inline complex pow(const complex & z, int power)
	{
		if (power < 0) {
			return pow(z, -power);
		}
		else if (power==0)
			return complex("1.0","0.0");
		else if(power==1)
			return z;
		else if(power==2)
			return z*z;
		else if(power==3)
			return z*z*z;
		else
		{
			unsigned int p = power;
			complex result("1.0","0.0"), z_to_the_current_power_of_two = z;
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
	
//	if(power==4)
//		return square(square(z));
//		else if(power==5)
//			return square(z)*cube(z);
//			else if(power==6)
//				return square(cube(z));
//				else
	
	
}




#endif




