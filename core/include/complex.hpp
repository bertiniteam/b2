#ifndef COMPLEX_CLASS_HPP
#define COMPLEX_CLASS_HPP



#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/random.hpp>


#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include <eigen3/Eigen/Core>

#include <assert.h>

// the following code block extends serialization to the mpfr_float class from boost::multiprecision
namespace boost { namespace serialization {
	
	
	/**
	 Save a mpfr_float type to a boost archive.
	 */
	template <typename Archive>
	void save(Archive& ar, ::boost::multiprecision::backends::mpfr_float_backend<0> const& r, unsigned /*version*/)
	{
		std::string tmp = r.str(0, std::ios::fixed);
		ar & tmp;
	}
	
	/**
	 Load a mpfr_float type from a boost archive.
	 */
	template <typename Archive>
	void load(Archive& ar, ::boost::multiprecision::backends::mpfr_float_backend<0>& r, unsigned /*version*/)
	{
		std::string tmp;
		ar & tmp;
		r = tmp.c_str();
	}
	
} } // re: namespaces

BOOST_SERIALIZATION_SPLIT_FREE(::boost::multiprecision::backends::mpfr_float_backend<0>)





namespace bertini {
	using boost::multiprecision::mpfr_float;
	
	class complex {
		
	private:
		// the real and imaginary parts of the complex number
		mpfr_float real_, imag_;
		
		
		// let the boost serialization library have access to the private members of this class.
		friend class boost::serialization::access;
		
		/**
		 \brief save method for archiving a bertini::complex
		 */
		template<class Archive>
		void save(Archive & ar, const unsigned int version) const
		{
			// note, version is always the latest when saving
			ar & real_;
			ar & imag_;
		}
		
		
		/**
		 \brief load method for archiving a bertini::complex
		 */
		template<class Archive>
		void load(Archive & ar, const unsigned int version)
		{
			ar & real_;
			ar & imag_;
		}
		
		BOOST_SERIALIZATION_SPLIT_MEMBER()
		// this has to be here so that the Boost.Serialization library
		// knows that there are different methods for the serialize() method.
		
		
	public:
		
		/**
		 default constructor
		 */
		complex():real_(), imag_(){}
		
		
		
		////////
		//
		//  construct real-valued complex numbers
		//
		//////////////
		
		/**
		 single-parameter for constructing a real-valued complex from a single real double number
		 */
		complex(double re) : real_(re), imag_("0.0"){}
		
		/**
		 single-parameter for constructing a real-valued complex from a single high-precision number
		 */
		complex(const mpfr_float & re) : real_(re), imag_("0.0"){}
		
		
		/**
		 single-parameter for constructing a real-valued complex from a convertible single string
		 */
		complex(const std::string & re) : real_(re), imag_("0.0"){}
		
		
		
		
		////////
		//
		//  construct complex numbers from two parameters
		//
		//////////////
		
		/**
		 two-parameter constructor for building a complex from two high precision numbers
		 */
		complex(const mpfr_float & re, const mpfr_float & im) : real_(re), imag_(im)
		{}
		
		
		/**
		 two-parameter constructor for building a complex from two low precision numbers
		 */
		complex(double re, double im) : real_(re), imag_(im)
		{}
		
		
		
		/**
		 two-parameter constructor for building a complex from two strings.
		 */
		complex(const std::string & re, const std::string & im) : real_(re), imag_(im)
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
		
		
		
		
		
		
		////////
		//
		//  other constructors
		//
		//////////////
		
		
//		/**
//		 the move constructor
//		 */
//		complex(complex&& other) : complex() // initialize via default constructor, C++11 only
//		{
//			swap(*this, other);
//		}
		
		
		
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
		
		
		
		
		
		
		
		////////
		//
		//  getters and setters
		//
		//////////////
		
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
		
		
		
		
		
		
		////////
		//
		//  some static functions which construct special numbers
		//
		//////////////
		
		
		
		
		
		inline static complex zero(){
			return complex("0.0","0.0");
		}
		
		inline static complex one(){
			return complex("1.0","0.0");
		}
		
		inline static complex two(){
			return complex("2.0","0.0");
		}
		
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
		
		
		
		
		
		
		
		
		
		
		////////
		//
		//  the fundamental arithmetic operators
		//
		//////////////
		
		
		
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
		 complex negation
		 */
		complex operator-() const
		{
			return complex(-real(), -imag());
		}
		
		
		
		
		
		
		
		
		
		
		
		/**
		 compute the square of the absolute value of the number
		 */
		mpfr_float abs2() const
		{
			return pow(real(),2)+pow(imag(),2);
		}
		
		/**
		 compute the absolute value of the number
		 */
		mpfr_float abs() const
		{
			return sqrt(abs2());
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
		 function may not tolerate white space in the number.
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
		
		
		
	}; // end declaration of the bertini::complex number class
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
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
		rhs.real(lhs - rhs.real());
		return rhs;
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
		return rhs*lhs; // it commutes!
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
	
	
	inline mpfr_float real(const complex & z)
	{
		return z.real();
	}
	
	inline mpfr_float imag(const complex & z)
	{
		return z.imag();
	}
	
	inline complex conj(const complex & z)
	{
		return z.conj();
	}
	
	
	
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
	 compute the inverse of a complex number
	 */
	inline complex inverse(const complex & z)
	{
		mpfr_float d = z.abs2();
		
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
	 
	 this could use fewer multiplications if it used more temporaries
	 */
	inline complex cube(const complex & z)
	{
		//		return complex(x^3 - 3*x*y^2, 3*x^2*y - y^3); // this deliberately left in for the equation.
		return complex(pow(z.real(),3) - mpfr_float("3.0")*z.real()*pow(z.imag(),2),
					   mpfr_float("3.0")*pow(z.real(),2)*z.imag() - pow(z.imag(),3));
		
	}
	
	
	
	/**
	 compute +,- integral powers of a complex number.
	 */
	inline complex pow(const complex & z, int power)
	{
		if (power < 0) {
			return pow(inverse(z), -power);
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
			unsigned int p(power);
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
	
	
	
	/**
	 construct a complex number from magnitude and angle.
	 */
	inline complex polar(const mpfr_float & rho, const mpfr_float & theta)
	{
		return complex(rho*cos(theta), rho*sin(theta));
	}
	
	
	/**
	 compute the square root of a complex number, using branch cut along the -x axis.
	 */
	inline complex sqrt(const complex & z)
	{
		return polar(sqrt(abs(z)), arg(z)/2);
	}
	
	
	/**
	 compute e^z for complex z.
	 */
	inline complex exp(const complex & z)
	{
		mpfr_float exp_of_x = exp(real(z));
		return complex(exp_of_x * cos(imag(z)), exp_of_x * sin(imag(z)));
	}
	
	
	inline complex sin(const complex & z)
	{
		return (exp(complex::i()*z) - exp(-complex::i()*z)) / complex::i() / mpfr_float("2.0");
	}
	
	inline complex cos(const complex & z)
	{
		return (exp(complex::i()*z) + exp(-complex::i()*z))  / mpfr_float("2.0");
	}
	
	inline complex tan(const complex & z)
	{
		return sin(z) / cos(z);
	}
	
	
	
	inline complex sinh(const complex & z)
	{
		return (exp(z) - exp(-z)) / mpfr_float("2.0");
	}
	
	inline complex cosh(const complex & z)
	{
		return (exp(z) + exp(-z))  / mpfr_float("2.0");
	}
	
	inline complex tanh(const complex & z)
	{
		return sinh(z) / cosh(z);
	}
	
	
	
	/**
	 complex logarithm base e.
	 */
	inline complex log(const complex & z)
	{
		return complex(log(abs(z)), arg(z));
	}
	
	
	
	/**
	 compute c^z, for c,z complex numbers
	 */
	inline complex pow(const complex & c, const complex & z)
	{
		return exp(c * log(z));
	}
	
	
	
	/**
	inverse sine of complex number
	 */
	inline complex asin(const complex & z)
	{
		return -complex::i() * log( complex::i()*z + sqrt( mpfr_float("1.0") - pow(z,2) ) );
	}
	
	
	/**
	 inverse cosine of complex number
	 */
	inline complex acos(const complex & z)
	{
		return -complex::i() * log( z + complex::i()*sqrt( mpfr_float("1.0") - pow(z,2) ) );
	}
	
	
	
	
	/**
	 inverse tangent of complex number
	 */
	inline complex atan(const complex & z)
	{
		return complex::i()/mpfr_float("2.0") * log( (complex::i() + z) / (complex::i() - z) );
	}
	
	
	
	
	/**
	 inverse hyperbolic sine of complex number
	 */
	inline complex asinh(const complex & z)
	{
		return log(
				   z + sqrt(
							square(z)+mpfr_float("1.0")
							)
				   );
	}
	
	/**
	 inverse hyperbolic cosine of complex number
	 */
	inline complex acosh(const complex & z)
	{
		return log(
				   z + sqrt(
							square(z)-mpfr_float("1.0")
							)
				   );
	}
	
	/**
	 inverse hyperbolic tangent of complex number
	 */
	inline complex atanh(const complex & z)
	{
		return mpfr_float("0.5") * log( (mpfr_float("1.0")+z)/(mpfr_float("1.0")-z) );
	}
	
	
	
	
} // re: namespace






// reopen the Eigen namespace to inject this struct.
namespace Eigen {
	
	using boost::multiprecision::mpfr_float;
	using namespace boost::multiprecision;
	template<> struct NumTraits<boost::multiprecision::mpfr_float> : GenericNumTraits<boost::multiprecision::mpfr_float> // permits to get the epsilon, dummy_precision, lowest, highest functions
	{
		
		typedef boost::multiprecision::mpfr_float Real;
		typedef boost::multiprecision::mpfr_float NonInteger;
		typedef boost::multiprecision::mpfr_float Nested;
		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned = 1,
			RequireInitialization = 1, // yes, require initialization, otherwise get crashes
			ReadCost = 20,
			AddCost = 30,
			MulCost = 40
		};
		

		inline static Real highest() {
			
			return (boost::multiprecision::mpfr_float(1) - epsilon()) * pow(boost::multiprecision::mpfr_float(2),mpfr_get_emax()-1);//);//boost::multiprecision::mpfr_float::default_precision());
			}
		
		inline static Real lowest() {
				return -highest();
			}
		
		inline static Real dummy_precision()
		{
			return pow( boost::multiprecision::mpfr_float(10),-int(boost::multiprecision::mpfr_float::default_precision()-3));
		}
		
		inline static Real epsilon()
		{
			return pow(boost::multiprecision::mpfr_float(10),-int(boost::multiprecision::mpfr_float::default_precision()));
		}
		//http://www.manpagez.com/info/mpfr/mpfr-2.3.2/mpfr_31.php
	};
	
	
	
	
	/**
	 \brief this templated struct permits us to use the Float type in Eigen matrices.
	 */
	template<> struct NumTraits<bertini::complex> : NumTraits<boost::multiprecision::mpfr_float> // permits to get the epsilon, dummy_precision, lowest, highest functions
	{
		typedef boost::multiprecision::mpfr_float Real;
		typedef boost::multiprecision::mpfr_float NonInteger;
		typedef bertini::complex Nested;// Nested;
		enum {
			IsComplex = 1,
			IsInteger = 0,
			IsSigned = 1,
			RequireInitialization = 1, // yes, require initialization, otherwise get crashes
			ReadCost = 2 * NumTraits<Real>::ReadCost,
			AddCost = 2 * NumTraits<Real>::AddCost,
			MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
		};
		
	};
	
	
	
	
	
	
	
}







#endif




