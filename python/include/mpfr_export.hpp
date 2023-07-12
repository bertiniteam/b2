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
//  silviana amethyst
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
#include "eigenpy_interaction.hpp"

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
			
			static T conj(T const& x){ return conj(x);}

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


	}
}



#endif
