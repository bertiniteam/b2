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
				.def("__eq__",&MPFRBaseVisitor::__eq__).def("__ne__",&MPFRBaseVisitor::__ne__)
				.def("__mul__",&MPFRBaseVisitor::__mul__).def("__imul__",&MPFRBaseVisitor::__imul__)
				.def("__div__",&MPFRBaseVisitor::__div__).def("__idiv__",&MPFRBaseVisitor::__idiv__)
				.def("__pow__",&MPFRBaseVisitor::__pow__)
				.def("__str__", &MPFRBaseVisitor::__str__)
				.def("__repr__", &MPFRBaseVisitor::__repr__)
				;
			}
			
			
		private:
			static bool __eq__(const MPFRBaseT& a, const MPFRBaseT& b)
			{
				return a == b;
			}
			
			static bool __ne__(const MPFRBaseT& a, const MPFRBaseT& b){ return !__eq__(a,b); }
			static MPFRBaseT __neg__(const MPFRBaseT& a){ return -a; };
			static MPFRBaseT __add__(const MPFRBaseT& a, const MPFRBaseT& b){ return a+b; };
			static MPFRBaseT __sub__(const MPFRBaseT& a, const MPFRBaseT& b){ return a-b; };
			static MPFRBaseT __iadd__(MPFRBaseT& a, const MPFRBaseT& b){ a+=b; return a; };
			static MPFRBaseT __isub__(MPFRBaseT& a, const MPFRBaseT& b){ a-=b; return a; };
			static MPFRBaseT __mul__(const MPFRBaseT& a, const MPFRBaseT& b){ return a*b; };
			static MPFRBaseT __div__(const MPFRBaseT& a, const MPFRBaseT& b){ return a/b; };
			static MPFRBaseT __imul__(MPFRBaseT& a, const MPFRBaseT& b){ a*=b; return a; };
			static MPFRBaseT __idiv__(MPFRBaseT& a, const MPFRBaseT& b){ a/=b; return a; };
			static MPFRBaseT __pow__(const MPFRBaseT& a, const MPFRBaseT& b){ return pow(a,b); };

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
