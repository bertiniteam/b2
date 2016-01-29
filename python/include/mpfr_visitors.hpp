//
//  mpfr_visitors.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 1/29/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_mpfr_visitors_hpp
#define Xcode_b2_mpfr_visitors_hpp
#include <boost/python.hpp>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
// #include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <boost/python/wrapper.hpp>

#include <boost/python/operators.hpp>
#include <boost/operators.hpp>


#include <bertini2/mpfr_extensions.hpp>

#include <sstream>



namespace bertini{
	namespace python{
		
		using namespace boost::python;
		
		template<typename MPFRBaseT>
		class MPFRBaseVisitor: public def_visitor<MPFRBaseVisitor<MPFRBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				cl
				.def(init<std::string>())
				.def(init<int>())
				.def(init<long int>())
				.def(init<MPFRBaseT>())
				
				.def("__neg__",&MPFRBaseVisitor::__neg__)
				
				.def("__add__",&MPFRBaseVisitor::__add__).def("__iadd__",&MPFRBaseVisitor::__iadd__)
				.def("__sub__",&MPFRBaseVisitor::__sub__).def("__isub__",&MPFRBaseVisitor::__isub__)
				.def("__mul__",&MPFRBaseVisitor::__mul__).def("__imul__",&MPFRBaseVisitor::__imul__)
				.def("__div__",&MPFRBaseVisitor::__div__).def("__idiv__",&MPFRBaseVisitor::__idiv__)
				.def("__pow__",&MPFRBaseVisitor::__pow__)

				.def("__add__",&MPFRBaseVisitor::__add_int).def("__iadd__",&MPFRBaseVisitor::__iadd_int)
				.def("__radd__",&MPFRBaseVisitor::__radd_int)
				.def("__sub__",&MPFRBaseVisitor::__sub_int).def("__isub__",&MPFRBaseVisitor::__isub_int)
				.def("__rsub__",&MPFRBaseVisitor::__rsub_int)
				.def("__mul__",&MPFRBaseVisitor::__mul_int).def("__imul__",&MPFRBaseVisitor::__imul_int)
				.def("__rmul__",&MPFRBaseVisitor::__rmul_int)
				.def("__div__",&MPFRBaseVisitor::__div_int).def("__idiv__",&MPFRBaseVisitor::__idiv_int)
				.def("__rdiv__",&MPFRBaseVisitor::__rdiv_int)
				.def("__pow__",&MPFRBaseVisitor::__pow_int)

				.def("__str__", &MPFRBaseVisitor::__str__)
				.def("__repr__", &MPFRBaseVisitor::__repr__)

				.def(self < self)
				.def(self <= self)
				.def(self > self)
				.def(self >= self)
				.def(self == self)
				.def(self != self)
				
				.def("precision", bmpprec1)
				.def("precision", bmpprec2)

				;
				
				def("default_precision", def_prec1);
				def("default_precision", def_prec2);
				
				
				def("abs", &MPFRBaseVisitor::__abs__);
				def("exp", &MPFRBaseVisitor::__exp__);
				def("log", &MPFRBaseVisitor::__log__);
				def("sqrt", &MPFRBaseVisitor::__sqrt__);

				def("sin", &MPFRBaseVisitor::__sin__);
				def("cos", &MPFRBaseVisitor::__cos__);
				def("tan", &MPFRBaseVisitor::__tan__);

				def("asin", &MPFRBaseVisitor::__asin__);
				def("acos", &MPFRBaseVisitor::__acos__);
				def("atan", &MPFRBaseVisitor::__atan__);

				def("sinh", &MPFRBaseVisitor::__sinh__);
				def("cosh", &MPFRBaseVisitor::__cosh__);
				def("tanh", &MPFRBaseVisitor::__tanh__);

			}
			
			
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
			}
			
			static std::string __repr__(const object& obj)
			{
				std::ostringstream oss;
				const MPFRBaseT& self=extract<MPFRBaseT>(obj)();
				std::stringstream ss;
				ss << self.str(0,std::ios::scientific);
				return ss.str();
			}

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

			unsigned (*def_prec1)() = &MPFRBaseT::default_precision;
			void (*def_prec2)(unsigned) = &MPFRBaseT::default_precision;
			
			unsigned (MPFRBaseT::*bmpprec1)() const= &MPFRBaseT::precision;
			void (MPFRBaseT::*bmpprec2)(unsigned) = &MPFRBaseT::precision;



		};
		
		
		
		template<typename MPFRBaseT>
		class MPFRFloatVisitor: public def_visitor<MPFRFloatVisitor<MPFRBaseT>>
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				MPFRBaseVisitor<MPFRBaseT>().visit(cl);
				
			}
			

			
			
		};
	}
}



#endif
