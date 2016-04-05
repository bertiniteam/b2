//This file is part of Bertini 2.0.
//
// python/bertini_python.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/bertini_python.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/bertini_python.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
//  python/mpfr_export.hpp:  Header file for exposing all multiprecision data types, those from boost and bertini::complex.

#ifndef Xcode_b2_mpfr_export_hpp
#define Xcode_b2_mpfr_export_hpp

#include "python_common.hpp"

namespace bertini{
	namespace python{
		
		using namespace boost::python;
		
		
		void ExportMpfr();
		
		
		// NOTE: Must redefine all aritmethic operators to support expresion templating in the mp data types
		
		
		
		
		/**
		 A base visitor exposing methods common to all multiprecision data types.  Mostly arithmetic and printing methods.
		 
		*/
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

			static MPFRBaseT __abs__(MPFRBaseT const& x){ return abs(x);}
			
			static std::string __str__(const object& obj)
			{
				std::ostringstream oss;
				const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
				std::stringstream ss;
				ss << self;
				return ss.str();
			};


		};

		
		
		
		
		/**
		 A base visitor exposing methods common to all float multiprecision data types.  These include trancendental functions and precision method.
		*/
		template<typename MPFRBaseT>
		class MPFRFloatBaseVisitor: public def_visitor<MPFRFloatBaseVisitor<MPFRBaseT> >
		{
			friend class def_visitor_access;

		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			static MPFRBaseT __div__(const MPFRBaseT& a, const MPFRBaseT& b){ return a/b; };
			static MPFRBaseT __idiv__(MPFRBaseT& a, const MPFRBaseT& b){ a/=b; return a; };
			static MPFRBaseT __pow__(const MPFRBaseT& a, const MPFRBaseT& b){ return pow(a,b); };
			
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


		
		
		
		/**
		 Visitor exposing important methods from the boost::multiprecision::mpz_int class
		*/
		template<typename MPFRBaseT>
		class MPFRIntVisitor: public def_visitor<MPFRIntVisitor<MPFRBaseT>>
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
			static MPFRBaseT __pow__(const MPFRBaseT& a, const int& b){ return pow(a,b); };
			
			static std::string __repr__(const object& obj)
			{
				std::ostringstream oss;
				const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
				std::stringstream ss;
				ss << self.str(0,std::ios::scientific);
				return ss.str();
			};
			
		};

	
		
		
		
		/**
		 Visitor exposing important methods from the boost::multiprecision::mpq_rational class
		 */
		template<typename MPFRBaseT>
		class MPFRRationalVisitor: public def_visitor<MPFRRationalVisitor<MPFRBaseT>>
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

			static MPFRBaseT __div__(const MPFRBaseT& a, const MPFRBaseT& b){ return a/b; };
			static MPFRBaseT __idiv__(MPFRBaseT& a, const MPFRBaseT& b){ a/=b; return a; };

			static MPFRBaseT __add_mpint(const MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ return a+b; };
			static MPFRBaseT __iadd_mpint(MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ a+=b; return a; };
			static MPFRBaseT __radd_mpint(const MPFRBaseT& b, const boost::multiprecision::mpz_int& a){ return a+b; };
			static MPFRBaseT __sub_mpint(const MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ return a-b; };
			static MPFRBaseT __isub_mpint(MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ a-=b; return a; };
			static MPFRBaseT __rsub_mpint(const MPFRBaseT& b, const boost::multiprecision::mpz_int& a){ return a-b; };
			static MPFRBaseT __mul_mpint(const MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ return a*b; };
			static MPFRBaseT __imul_mpint(MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ a*=b; return a; };
			static MPFRBaseT __rmul_mpint(const MPFRBaseT& b, const boost::multiprecision::mpz_int& a){ return a*b; };
			static MPFRBaseT __div_mpint(const MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ return a/b; };
			static MPFRBaseT __idiv_mpint(MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ a/=b; return a; };

			static std::string __repr__(const object& obj)
			{
				std::ostringstream oss;
				const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
				std::stringstream ss;
				ss << self.str(0,std::ios::scientific);
				return ss.str();
			};
			
		};

		
		
		
		
		
		/**
		 Visitor exposing important methods from the bmp class, where bmp is defined in the bertini core.  This is the boost multiprecision real float class.
		 */
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
		
		
		
		
		
		
		
		
		/**
		 Visitor exposing important methods from the bertini::complex
		 */
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
