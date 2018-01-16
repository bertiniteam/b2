//This file is part of Bertini 2.
//
//python/mpfr_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/mpfr_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/mpfr_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.
//
// individual authors of this file include:
//
//  Danielle Brake
//  University of Wisconsin - Eau Claire
//  Fall 2017, Spring 2018
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
//  python/mpfr_export.hpp:  Header file for exposing all multiprecision data types, those from boost and bertini::complex.

#pragma once
#ifndef BERTINI_PYTHON_MPFR_EXPORT_HPP
#define BERTINI_PYTHON_MPFR_EXPORT_HPP

#include "python_common.hpp"

namespace bertini{
	namespace python{
		
		using namespace boost::python;
		
		
		void ExportMpfr();
		
		/**
		 \brief Exposes  precision
		*/
		template<typename T>
		class PrecisionVisitor: public def_visitor<PrecisionVisitor<T>>
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		private:

			unsigned (T::*get_prec)() const= &T::precision; // return type needs to be PrecT
			void (T::*set_prec)(unsigned) = &T::precision;
		};



		/**
		 \brief Exposes str, repr, and precision
		*/
		template<typename T>
		class RealStrVisitor: public def_visitor<RealStrVisitor<T>>
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		private:

			static std::string __str__(const object& obj)
			{
				std::ostringstream oss;
				const T& self=extract<T>(obj)();
				std::stringstream ss;
				ss << self;
				return ss.str();
			};

			static std::string __repr__(const object& obj)
			{
				std::ostringstream oss;
				const T& self=extract<T>(obj)();
				std::stringstream ss;
				ss << self.str(0,std::ios::scientific);
				return ss.str();
			};
		};



		/**
		 \brief Exposes == and != for homogeneous comparison
		*/
		template<typename T>
		class EqualitySelfVisitor: public def_visitor<EqualitySelfVisitor<T>>
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		private:
		};

		/**
		 \brief Exposes == and != for inhomogeneous comparison
		*/
		template<typename T, typename S>
		class EqualityVisitor: public def_visitor<EqualityVisitor<T,S>>
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		private:
		};

		/**
		 \brief Exposes +, - and *
		*/
		template<typename T, typename S>
		class RingVisitor: public def_visitor<RingVisitor<T, S> >
		{
			static_assert(!std::is_same<T,S>::value, "RingVisitor is to define T-S operations.  for T-T operations, use RingSelfVisitor");
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			static T __add__(const T& a, const S& b){ return a+b; };
			static T __radd__(const T& b, const S& a){ return a+b; };
			static T __iadd__(T& a, const S& b){ a+=b; return a; };

			static T __sub__(const T& a, const S& b){ return a-b; };
			static T __rsub__(const T& b, const S& a){ return a-b; };
			static T __isub__(T& a, const S& b){ a-=b; return a; };

			static T __mul__(const T& a, const S& b){ return a*b; };
			static T __rmul__(const T& b, const S& a){ return a*b; };
			static T __imul__(T& a, const S& b){ a*=b; return a; };
		}; // RingVisitor
		
		/**
		 \brief Exposes +, - and *
		*/
		template<typename T>
		class RingSelfVisitor: public def_visitor<RingSelfVisitor<T> >
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			static T __add__(const T& a, const T& b){ return a+b; };
			static T __iadd__(T& a, const T& b){ a+=b; return a; };

			static T __sub__(const T& a, const T& b){ return a-b; };
			static T __isub__(T& a, const T& b){ a-=b; return a; };

			static T __mul__(const T& a, const T& b){ return a*b; };
			static T __imul__(T& a, const T& b){ a*=b; return a; };

			static T __neg__(const T& a){ return -a; };
		}; // RingSelfVisitor
		
		/**
		 \brief Exposes +, - and *
		*/
		template<typename T>
		class RealFreeVisitor: public def_visitor<RealFreeVisitor<T> >
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			static T __abs__(T const& x){ return abs(x);}
		}; // RingSelfVisitor
		


		/**
		 \brief Exposes +,-,*,/
		*/
		template<typename T, typename S>
		class FieldVisitor: public def_visitor<FieldVisitor<T, S> >
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			static T __div__(const T& a, const S& b){ return a/b; };
			static T __rdiv__(const T& b, const S& a){ return a/b; };
			static T __idiv__(T& a, const S& b){ a/=b; return a; };
		}; // FieldVisitor

		template<typename T>
		class FieldSelfVisitor: public def_visitor<FieldSelfVisitor<T> >
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			static T div(const T& a, const T& b){ return a/b; };
			static T idiv(T& a, const T& b){ a/=b; return a; };
		}; // FieldSelfVisitor



		template<typename T, typename S>
		class PowVisitor: public def_visitor<PowVisitor<T,S>>
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:
			static T __pow__(const T& a, const S& b){using std::pow; using boost::multiprecision::pow; return pow(a,b); };
		};


		template<typename T, typename S>
		class GreatLessVisitor: public def_visitor<GreatLessVisitor<T,S>>
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:

		};

			
		template<typename T>
		class GreatLessSelfVisitor: public def_visitor<GreatLessSelfVisitor<T>>
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:

		};



		template<typename T>
		class TranscendentalVisitor: public def_visitor<TranscendentalVisitor<T>>
		{
			friend class def_visitor_access;
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;

		private:
		
			static T __log__(T const& x){ return log(x);}
			static T __exp__(T const& x){ return exp(x);}
			static T __sqrt__(T const& x){ return sqrt(x);}
			
			static T __sin__(T const& x){ return sin(x);}
			static T __cos__(T const& x){ return cos(x);}
			static T __tan__(T const& x){ return tan(x);}
			
			static T __asin__(T const& x){ return asin(x);}
			static T __acos__(T const& x){ return acos(x);}
			static T __atan__(T const& x){ return atan(x);}
			
			static T __sinh__(T const& x){ return sinh(x);}
			static T __cosh__(T const& x){ return cosh(x);}
			static T __tanh__(T const& x){ return tanh(x);}

			static T __asinh__(T const& x){ return asinh(x);}
			static T __acosh__(T const& x){ return acosh(x);}
			static T __atanh__(T const& x){ return atanh(x);}
		};

		/**
		 \brief Exposes complex-specific things
		 */
		template<typename T>
		class ComplexVisitor: public def_visitor<ComplexVisitor<T>>
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			using RealT = typename NumTraits<T>::Real;
			
		private:
			static void set_real(T &c, mpfr_float const& r) { c.real(r);}
			static RealT get_real(T const&c) { return c.real();}
			
			static void set_imag(T &c, RealT const& r) { c.imag(r);}
			static RealT get_imag(T const&c) { return c.imag();}

			static RealT __abs__(T const& x){ return abs(x);}
			
			static std::string __str__(const object& obj)
			{
				std::ostringstream oss;
				const T& self=extract<T>(obj)();
				std::stringstream ss;
				ss << self;
				return ss.str();
			}
			
			static std::string __repr__(const object& obj)
			{
				std::ostringstream oss;
				const T& self=extract<T>(obj)();
				std::stringstream ss;
				ss << "(" << real(self).str(0,std::ios::scientific) << ", " << imag(self).str(0,std::ios::scientific) << ")";
				return ss.str();
			}

		};
		// /**
		//  A base visitor exposing methods common to all multiprecision data types.  Mostly arithmetic and printing methods.
		 
		// */
		// template<typename MPFRBaseT>
		// class MPFRBaseVisitor: public def_visitor<MPFRBaseVisitor<MPFRBaseT> >
		// {
		// public:
		// 	template<class PyClass>
		// 	void visit(PyClass& cl) const;
			
			
		// private:
		// 	static MPFRBaseT __neg__(const MPFRBaseT& a){ return -a; };
			
		// 	static MPFRBaseT __add__(const MPFRBaseT& a, const MPFRBaseT& b){ return a+b; };
		// 	static MPFRBaseT __iadd__(MPFRBaseT& a, const MPFRBaseT& b){ a+=b; return a; };
		// 	static MPFRBaseT __sub__(const MPFRBaseT& a, const MPFRBaseT& b){ return a-b; };
		// 	static MPFRBaseT __isub__(MPFRBaseT& a, const MPFRBaseT& b){ a-=b; return a; };
		// 	static MPFRBaseT __mul__(const MPFRBaseT& a, const MPFRBaseT& b){ return a*b; };
		// 	static MPFRBaseT __imul__(MPFRBaseT& a, const MPFRBaseT& b){ a*=b; return a; };

		// 	static MPFRBaseT __abs__(MPFRBaseT const& x){ return abs(x);}
			
		// 	static std::string __str__(const object& obj)
		// 	{
		// 		std::ostringstream oss;
		// 		const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
		// 		std::stringstream ss;
		// 		ss << self;
		// 		return ss.str();
		// 	};


		// };

		
		
		
		
		// /**
		//  A base visitor exposing methods common to all float multiprecision data types.  These include trancendental functions and precision method.
		// */
		// template<typename MPFRBaseT>
		// class MPFRFloatBaseVisitor: public def_visitor<MPFRFloatBaseVisitor<MPFRBaseT> >
		// {
		// 	friend class def_visitor_access;

		// public:
		// 	template<class PyClass>
		// 	void visit(PyClass& cl) const;
			
			
		// private:
		// 	static MPFRBaseT __div__(const MPFRBaseT& a, const MPFRBaseT& b){ return a/b; };
		// 	static MPFRBaseT __idiv__(MPFRBaseT& a, const MPFRBaseT& b){ a/=b; return a; };
		// 	static MPFRBaseT __pow__(const MPFRBaseT& a, const MPFRBaseT& b){ return pow(a,b); };
			
		// 	static MPFRBaseT __log__(MPFRBaseT const& x){ return log(x);}
		// 	static MPFRBaseT __exp__(MPFRBaseT const& x){ return exp(x);}
		// 	static MPFRBaseT __sqrt__(MPFRBaseT const& x){ return sqrt(x);}
			
		// 	static MPFRBaseT __sin__(MPFRBaseT const& x){ return sin(x);}
		// 	static MPFRBaseT __cos__(MPFRBaseT const& x){ return cos(x);}
		// 	static MPFRBaseT __tan__(MPFRBaseT const& x){ return tan(x);}
			
		// 	static MPFRBaseT __asin__(MPFRBaseT const& x){ return asin(x);}
		// 	static MPFRBaseT __acos__(MPFRBaseT const& x){ return acos(x);}
		// 	static MPFRBaseT __atan__(MPFRBaseT const& x){ return atan(x);}
			
		// 	static MPFRBaseT __sinh__(MPFRBaseT const& x){ return sinh(x);}
		// 	static MPFRBaseT __cosh__(MPFRBaseT const& x){ return cosh(x);}
		// 	static MPFRBaseT __tanh__(MPFRBaseT const& x){ return tanh(x);}
			
		// 	unsigned (MPFRBaseT::*bmpprec1)() const= &MPFRBaseT::precision;
		// 	void (MPFRBaseT::*bmpprec2)(unsigned) = &MPFRBaseT::precision;
			
			
			
		// };


		
		
		
		// /**
		//  Visitor exposing important methods from the boost::multiprecision::mpz_int class
		// */
		// template<typename T>
		// class IntVisitor: public def_visitor<IntVisitor<T>>
		// {
		// 	friend class def_visitor_access;
			
		// public:
		// 	template<class PyClass>
		// 	void visit(PyClass& cl) const;
			
		// private:
			
		// };

	
		
		
		
		// /**
		//  Visitor exposing important methods from the boost::multiprecision::mpq_rational class
		//  */
		// template<typename MPFRBaseT>
		// class MPFRRationalVisitor: public def_visitor<MPFRRationalVisitor<MPFRBaseT>>
		// {
		// 	friend class def_visitor_access;
			
		// public:
		// 	template<class PyClass>
		// 	void visit(PyClass& cl) const;
			
			
			
		// private:
		// 	static MPFRBaseT __add_int(const MPFRBaseT& a, const int& b){ return a+b; };
		// 	static MPFRBaseT __iadd_int(MPFRBaseT& a, const int& b){ a+=b; return a; };
		// 	static MPFRBaseT __radd_int(const MPFRBaseT& b, const int& a){ return a+b; };
		// 	static MPFRBaseT __sub_int(const MPFRBaseT& a, const int& b){ return a-b; };
		// 	static MPFRBaseT __isub_int(MPFRBaseT& a, const int& b){ a-=b; return a; };
		// 	static MPFRBaseT __rsub_int(const MPFRBaseT& b, const int& a){ return a-b; };
		// 	static MPFRBaseT __mul_int(const MPFRBaseT& a, const int& b){ return a*b; };
		// 	static MPFRBaseT __imul_int(MPFRBaseT& a, const int& b){ a*=b; return a; };
		// 	static MPFRBaseT __rmul_int(const MPFRBaseT& b, const int& a){ return a*b; };
		// 	static MPFRBaseT __div_int(const MPFRBaseT& a, const int& b){ return a/b; };
		// 	static MPFRBaseT __idiv_int(MPFRBaseT& a, const int& b){ a/=b; return a; };
		// 	static MPFRBaseT __rdiv_int(const MPFRBaseT& b, const int& a){ return a/b; };

		// 	static MPFRBaseT __div__(const MPFRBaseT& a, const MPFRBaseT& b){ return a/b; };
		// 	static MPFRBaseT __idiv__(MPFRBaseT& a, const MPFRBaseT& b){ a/=b; return a; };

		// 	static MPFRBaseT __add_mpint(const MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ return a+b; };
		// 	static MPFRBaseT __iadd_mpint(MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ a+=b; return a; };
		// 	static MPFRBaseT __radd_mpint(const MPFRBaseT& b, const boost::multiprecision::mpz_int& a){ return a+b; };
		// 	static MPFRBaseT __sub_mpint(const MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ return a-b; };
		// 	static MPFRBaseT __isub_mpint(MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ a-=b; return a; };
		// 	static MPFRBaseT __rsub_mpint(const MPFRBaseT& b, const boost::multiprecision::mpz_int& a){ return a-b; };
		// 	static MPFRBaseT __mul_mpint(const MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ return a*b; };
		// 	static MPFRBaseT __imul_mpint(MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ a*=b; return a; };
		// 	static MPFRBaseT __rmul_mpint(const MPFRBaseT& b, const boost::multiprecision::mpz_int& a){ return a*b; };
		// 	static MPFRBaseT __div_mpint(const MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ return a/b; };
		// 	static MPFRBaseT __idiv_mpint(MPFRBaseT& a, const boost::multiprecision::mpz_int& b){ a/=b; return a; };

		// 	static std::string __repr__(const object& obj)
		// 	{
		// 		std::ostringstream oss;
		// 		const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
		// 		std::stringstream ss;
		// 		ss << self.str(0,std::ios::scientific);
		// 		return ss.str();
		// 	};
			
		// };

		
		
		
		
		
		// /**
		//  Visitor exposing important methods from the bmp class, where bmp is defined in the bertini core.  This is the boost multiprecision real float class.
		//  */
		// template<typename MPFRBaseT>
		// class MPFRFloatVisitor: public def_visitor<MPFRFloatVisitor<MPFRBaseT>>
		// {
		// 	friend class def_visitor_access;
			
		// public:
		// 	template<class PyClass>
		// 	void visit(PyClass& cl) const;
			
	
			
		// private:
		// 	static MPFRBaseT __add_int(const MPFRBaseT& a, const int& b){ return a+b; };
		// 	static MPFRBaseT __iadd_int(MPFRBaseT& a, const int& b){ a+=b; return a; };
		// 	static MPFRBaseT __radd_int(const MPFRBaseT& b, const int& a){ return a+b; };
		// 	static MPFRBaseT __sub_int(const MPFRBaseT& a, const int& b){ return a-b; };
		// 	static MPFRBaseT __isub_int(MPFRBaseT& a, const int& b){ a-=b; return a; };
		// 	static MPFRBaseT __rsub_int(const MPFRBaseT& b, const int& a){ return a-b; };
		// 	static MPFRBaseT __mul_int(const MPFRBaseT& a, const int& b){ return a*b; };
		// 	static MPFRBaseT __imul_int(MPFRBaseT& a, const int& b){ a*=b; return a; };
		// 	static MPFRBaseT __rmul_int(const MPFRBaseT& b, const int& a){ return a*b; };
		// 	static MPFRBaseT __div_int(const MPFRBaseT& a, const int& b){ return a/b; };
		// 	static MPFRBaseT __idiv_int(MPFRBaseT& a, const int& b){ a/=b; return a; };
		// 	static MPFRBaseT __rdiv_int(const MPFRBaseT& b, const int& a){ return a/b; };
		// 	static MPFRBaseT __pow_int(const MPFRBaseT& a, const int& b){ return pow(a,b); };
			
			
		// 	static std::string __repr__(const object& obj)
		// 	{
		// 		std::ostringstream oss;
		// 		const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
		// 		std::stringstream ss;
		// 		ss << self.str(0,std::ios::scientific);
		// 		return ss.str();
		// 	};
			
		// 	unsigned (*def_prec1)() = &MPFRBaseT::default_precision;
		// 	void (*def_prec2)(unsigned) = &MPFRBaseT::default_precision;
		// };
		
		
		
		
		
		
		
		
		// /**
		//  Visitor exposing important methods from the bertini::complex
		//  */
		// template<typename MPFRBaseT>
		// class MPFRComplexVisitor: public def_visitor<MPFRComplexVisitor<MPFRBaseT>>
		// {
		// 	friend class def_visitor_access;
			
		// public:
		// 	template<class PyClass>
		// 	void visit(PyClass& cl) const;
			
			
		// private:
		// 	template<typename T>
		// 	static MPFRBaseT __add_float(const MPFRBaseT& a, const T& b){ return a+b; };
		// 	static MPFRBaseT __iadd_float(MPFRBaseT& a, const mpfr_float& b){ a+=b; return a; };
		// 	static MPFRBaseT __radd_float(const MPFRBaseT& b, const mpfr_float& a){ return a+b; };
		// 	static MPFRBaseT __sub_float(const MPFRBaseT& a, const mpfr_float& b){ return a-b; };
		// 	static MPFRBaseT __isub_float(MPFRBaseT& a, const mpfr_float& b){ a-=b; return a; };
		// 	static MPFRBaseT __rsub_float(const MPFRBaseT& b, const mpfr_float& a){ return a-b; };
		// 	static MPFRBaseT __mul_float(const MPFRBaseT& a, const mpfr_float& b){ return a*b; };
		// 	static MPFRBaseT __imul_float(MPFRBaseT& a, const mpfr_float& b){ a*=b; return a; };
		// 	static MPFRBaseT __rmul_float(const MPFRBaseT& b, const mpfr_float& a){ return a*b; };
		// 	static MPFRBaseT __div_float(const MPFRBaseT& a, const mpfr_float& b){ return a/b; };
		// 	static MPFRBaseT __idiv_float(MPFRBaseT& a, const mpfr_float& b){ a/=b; return a; };
		// 	static MPFRBaseT __rdiv_float(const MPFRBaseT& b, const mpfr_float& a){ return a/b; };
		// 	static MPFRBaseT __pow_int(const MPFRBaseT& a, const int& b){ return pow(a,b); };
		// 	static MPFRBaseT __pow_float(const MPFRBaseT& a, const mpfr_float& b){ return pow(a,b); };
			
		// 	static std::string __str__(const object& obj)
		// 	{
		// 		std::ostringstream oss;
		// 		const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
		// 		std::stringstream ss;
		// 		ss << self;
		// 		return ss.str();
		// 	}
			
		// 	static std::string __repr__(const object& obj)
		// 	{
		// 		std::ostringstream oss;
		// 		const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
		// 		std::stringstream ss;
		// 		ss << "(" << real(self).str(0,std::ios::scientific) << ", " << imag(self).str(0,std::ios::scientific) << ")";
		// 		return ss.str();
		// 	}

		// 	static void set_real(MPFRBaseT &c, mpfr_float const& r) { c.real(r);}
		// 	static mpfr_float get_real(MPFRBaseT const&c) { return c.real();}
			
		// 	static void set_imag(MPFRBaseT &c, mpfr_float const& r) { c.imag(r);}
		// 	static mpfr_float get_imag(MPFRBaseT const&c) { return c.imag();}
		// };

	}
}



#endif
