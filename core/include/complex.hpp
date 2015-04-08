#ifndef COMPLEX_CLASS_HPP
#define COMPLEX_CLASS_HPP

#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>

#include <assert.h>

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
		
		
		inline mpfr_float real() const {return real_;}
		inline mpfr_float imag() const {return imag_;}
		
		void real(const mpfr_float & new_real){real_ = new_real;}
		void imag(const mpfr_float & new_imag){imag_ = new_imag;}
		
		void real(const std::string & new_real){real_ = mpfr_float(new_real);}
		void imag(const std::string & new_imag){imag_ = mpfr_float(new_imag);}
		
		
		complex& operator+=(const complex & rhs)
		{
			real_+=rhs.real_;
			imag_+=rhs.imag_;
			return *this;
		}
		
		
		
		complex& operator-=(const complex & rhs)
		{
			real_-=rhs.real_;
			imag_-=rhs.imag_;
			return *this;
		}
		
		
		/**
		 complex multiplication.  uses a single temporary variable
		 */
		complex& operator*=(const complex & rhs)
		{
			mpfr_float a = real_*rhs.real_ - imag_*rhs.imag_;
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
			mpfr_float a = real_*rhs.real_ + imag_*rhs.imag_;
			imag_ = imag_*rhs.real_ - real_*rhs.imag_/d;
			real_ = a/d;
			
			return *this;
		}
		
		
		
		
		
		
		
		
		
		friend std::ostream& operator<<(std::ostream& out, const complex & z)
		{
			out << "(" << z.real() << "," << z.imag() << ")";
			return out;
		}
		
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
		
		
		
		
		mpfr_float arg() const
		{
			return boost::multiprecision::atan2(imag(),real());
		}
		
		
		
		mpfr_float norm() const
		{
			return abs2();
		}
		
		
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
	
	
	inline mpfr_float abs2(const complex & z)
	{
		return z.abs2();
	}
	
	inline mpfr_float abs(const complex & z)
	{
		return boost::multiprecision::sqrt(abs2(z));
	}
	
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
	
	
}


#endif


