//
//  mpfr_visitors.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 1/29/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_mpfr_visitors_hpp
#define Xcode_b2_mpfr_visitors_hpp

#include "python_common.hpp"

namespace bertini{
	namespace python{
		
		using namespace boost::python;
		//typedef boost::multiprecision::number<boost::multiprecision::mpfr_float, boost::multiprecision::et_off> bmp;
		//		typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0> , boost::multiprecision::et_off> bmp;
		
		
		void ExportMpfr();
		
		template<typename MPFRBaseT>
		class MPFRBaseVisitor: public def_visitor<MPFRBaseVisitor<MPFRBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			static MPFRBaseT __neg__(const MPFRBaseT& a){ return -a; };
			
			static MPFRBaseT __add__(const MPFRBaseT& a, const MPFRBaseT& b){ return a+b; };
			static MPFRBaseT __iadd__(MPFRBaseT& a, const MPFRBaseT& b){ a+=b; return a; };
			static MPFRBaseT __sub__(const MPFRBaseT& a, const MPFRBaseT& b){ return a-b; };
			static MPFRBaseT __isub__(MPFRBaseT& a, const MPFRBaseT& b){ a-=b; return a; };
			static MPFRBaseT __mul__(const MPFRBaseT& a, const MPFRBaseT& b){ return a*b; };
			static MPFRBaseT __imul__(MPFRBaseT& a, const MPFRBaseT& b){ a*=b; return a; };
			static MPFRBaseT __div__(const MPFRBaseT& a, const MPFRBaseT& b){ return a/b; };
			static MPFRBaseT __idiv__(MPFRBaseT& a, const MPFRBaseT& b){ a/=b; return a; };
			static MPFRBaseT __pow__(const MPFRBaseT& a, const MPFRBaseT& b){ return pow(a,b); };

			static MPFRBaseT __abs__(MPFRBaseT const& x){ return abs(x);}
			static MPFRBaseT __log__(MPFRBaseT const& x){ return log(x);}
			static MPFRBaseT __exp__(MPFRBaseT const& x){ return exp(x);}
			static MPFRBaseT __sqrt__(MPFRBaseT const& x){ return sqrt(x);}

			static MPFRBaseT __sin__(MPFRBaseT const& x){ return sin(x);}
			static MPFRBaseT __cos__(MPFRBaseT const& x){ return cos(x);}
			static MPFRBaseT __tan__(MPFRBaseT const& x){ return tan(x);}

			static MPFRBaseT __asin__(MPFRBaseT const& x){ return asin(x);}
			static MPFRBaseT __acos__(MPFRBaseT const& x){ return acos(x);}
			static MPFRBaseT __atan__(MPFRBaseT const& x){ return atan(x);}

			static MPFRBaseT __sinh__(MPFRBaseT const& x){ return sinh(x);}
			static MPFRBaseT __cosh__(MPFRBaseT const& x){ return cosh(x);}
			static MPFRBaseT __tanh__(MPFRBaseT const& x){ return tanh(x);}
			
			unsigned (MPFRBaseT::*bmpprec1)() const= &MPFRBaseT::precision;
			void (MPFRBaseT::*bmpprec2)(unsigned) = &MPFRBaseT::precision;



		};
		
		
		// boost::multiprecision float class
		template<typename MPFRBaseT>
		class MPFRFloatVisitor: public def_visitor<MPFRFloatVisitor<MPFRBaseT>>
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
	
			
		private:
			static MPFRBaseT __add_int(const MPFRBaseT& a, const int& b){ return a+b; };
			static MPFRBaseT __iadd_int(MPFRBaseT& a, const int& b){ a+=b; return a; };
			static MPFRBaseT __radd_int(const MPFRBaseT& b, const int& a){ return a+b; };
			static MPFRBaseT __sub_int(const MPFRBaseT& a, const int& b){ return a-b; };
			static MPFRBaseT __isub_int(MPFRBaseT& a, const int& b){ a-=b; return a; };
			static MPFRBaseT __rsub_int(const MPFRBaseT& b, const int& a){ return a-b; };
			static MPFRBaseT __mul_int(const MPFRBaseT& a, const int& b){ return a*b; };
			static MPFRBaseT __imul_int(MPFRBaseT& a, const int& b){ a*=b; return a; };
			static MPFRBaseT __rmul_int(const MPFRBaseT& b, const int& a){ return a*b; };
			static MPFRBaseT __div_int(const MPFRBaseT& a, const int& b){ return a/b; };
			static MPFRBaseT __idiv_int(MPFRBaseT& a, const int& b){ a/=b; return a; };
			static MPFRBaseT __rdiv_int(const MPFRBaseT& b, const int& a){ return a/b; };
			static MPFRBaseT __pow_int(const MPFRBaseT& a, const int& b){ return pow(a,b); };
			
			static std::string __str__(const object& obj)
			{
				std::ostringstream oss;
				const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
				std::stringstream ss;
				ss << self;
				return ss.str();
			};
			
			static std::string __repr__(const object& obj)
			{
				std::ostringstream oss;
				const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
				std::stringstream ss;
				ss << self.str(0,std::ios::scientific);
				return ss.str();
			};
			
			unsigned (*def_prec1)() = &MPFRBaseT::default_precision;
			void (*def_prec2)(unsigned) = &MPFRBaseT::default_precision;
		};
		
		
		
		
		// bertini::complex class
		template<typename MPFRBaseT>
		class MPFRComplexVisitor: public def_visitor<MPFRComplexVisitor<MPFRBaseT>>
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			static MPFRBaseT __add_float(const MPFRBaseT& a, const mpfr_float& b){ return a+b; };
			static MPFRBaseT __iadd_float(MPFRBaseT& a, const mpfr_float& b){ a+=b; return a; };
			static MPFRBaseT __radd_float(const MPFRBaseT& b, const mpfr_float& a){ return a+b; };
			static MPFRBaseT __sub_float(const MPFRBaseT& a, const mpfr_float& b){ return a-b; };
			static MPFRBaseT __isub_float(MPFRBaseT& a, const mpfr_float& b){ a-=b; return a; };
			static MPFRBaseT __rsub_float(const MPFRBaseT& b, const mpfr_float& a){ return a-b; };
			static MPFRBaseT __mul_float(const MPFRBaseT& a, const mpfr_float& b){ return a*b; };
			static MPFRBaseT __imul_float(MPFRBaseT& a, const mpfr_float& b){ a*=b; return a; };
			static MPFRBaseT __rmul_float(const MPFRBaseT& b, const mpfr_float& a){ return a*b; };
			static MPFRBaseT __div_float(const MPFRBaseT& a, const mpfr_float& b){ return a/b; };
			static MPFRBaseT __idiv_float(MPFRBaseT& a, const mpfr_float& b){ a/=b; return a; };
			static MPFRBaseT __rdiv_float(const MPFRBaseT& b, const mpfr_float& a){ return a/b; };
			static MPFRBaseT __pow_int(const MPFRBaseT& a, const int& b){ return pow(a,b); };
			static MPFRBaseT __pow_float(const MPFRBaseT& a, const mpfr_float& b){ return pow(a,b); };
			
			static std::string __str__(const object& obj)
			{
				std::ostringstream oss;
				const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
				std::stringstream ss;
				ss << self;
				return ss.str();
			}
			
			static std::string __repr__(const object& obj)
			{
				std::ostringstream oss;
				const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
				std::stringstream ss;
				ss << "(" << real(self).str(0,std::ios::scientific) << ", " << imag(self).str(0,std::ios::scientific) << ")";
				return ss.str();
			}

			mpfr_float (MPFRBaseT::*getreal)() const = &MPFRBaseT::real;
			void (MPFRBaseT::*setreal)(const mpfr_float&) = &MPFRBaseT::real;
			mpfr_float (MPFRBaseT::*getimag)() const = &MPFRBaseT::imag;
			void (MPFRBaseT::*setimag)(const mpfr_float&) = &MPFRBaseT::imag;
			
			void set_real(MPFRBaseT &c, mpfr_float const& r) { c.real(r);}
			mpfr_float get_real(MPFRBaseT const&c) { return c.real();}
			
			void set_imag(MPFRBaseT &c, mpfr_float const& r) { c.imag(r);}
			mpfr_float get_imag(MPFRBaseT const&c) { return c.imag();}
		};

	}
}



#endif
