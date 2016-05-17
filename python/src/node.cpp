//This file is part of Bertini 2.
//
//python/node.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/node.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/node.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
//
//  python/node.cpp:  the source file for the python interface for Node class.

#include <stdio.h>

#include "node.hpp"


using namespace boost::python;

using dbl = std::complex<double>;
using mpfr = bertini::complex;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NodeEvalOverloadsDbl, bertini::node::Node::template Eval<dbl>, 0, 1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NodeEvalOverloadsMPFR, bertini::node::Node::template Eval<mpfr>, 0, 1);


namespace bertini{
	namespace python{
		void ExportNode()
		{
			using namespace bertini::node;
			using NodePtr = std::shared_ptr<node::Node>;
			

			
			class_<bertini::python::NodeWrap, boost::noncopyable, std::shared_ptr<node::Node> >("Node", no_init)
			.def("reset", &bertini::node::Node::Reset, &bertini::python::NodeWrap::default_Reset)
			
			.def("evald", &bertini::node::Node::Eval<dbl>,
				 NodeEvalOverloadsDbl(args("DiffVar"), "Eval's docstring"))
			.def("evalmp", &bertini::node::Node::Eval<mpfr>,
				 NodeEvalOverloadsMPFR(args("DiffVar"), "Eval's docstring"))
			.def(self_ns::str(self_ns::self))
			.def(self_ns::repr(self_ns::self))
			
			.def("__add__",addNodeNode)
			.def("__add__",addNodeDouble)
			.def("__radd__",addNodeDouble)
			.def("__add__",addNodeDbl)
			.def("__radd__",addNodeDbl)
			.def("__add__",addNodeMpfr)
			.def("__radd__",addNodeMpfr)
			.def("__add__",addNodeInt)
			.def("__radd__",addNodeInt)
			.def("__iadd__",&iaddNodeNode)
			.def("__iadd__",&iaddNodeDouble)
			.def("__iadd__", &iaddSumNode)

			.def("__sub__",subNodeNode)
			.def("__sub__",subNodeDouble)
			.def("__rsub__",subNodeDouble)
			.def("__sub__",subNodeDbl)
			.def("__rsub__",subNodeDbl)
			.def("__sub__",subNodeMpfr)
			.def("__rsub__",subNodeMpfr)
			.def("__sub__",subNodeInt)
			.def("__rsub__",subNodeInt)
			.def("__isub__",&isubNodeNode)
			.def("__isub__",&isubNodeDouble)
			.def("__isub__", &isubSumNode)

			.def("__mul__",multNodeNode)
			.def("__mul__",multNodeDouble)
			.def("__rmul__",multNodeDouble)
			.def("__mul__",multNodeDbl)
			.def("__rmul__",multNodeDbl)
			.def("__mul__",multNodeMpfr)
			.def("__rmul__",multNodeMpfr)
			.def("__mul__",multNodeInt)
			.def("__rmul__",multNodeInt)
			.def("__imul__",&imultNodeNode)
			.def("__imul__",&imultNodeDouble)
			.def("__imul__",imultMultNode)

			.def("__div__",divNodeNode)
			.def("__div__",divNodeDouble)
			.def("__rdiv__",divNodeDouble)
			.def("__div__",divNodeDbl)
			.def("__rdiv__",divNodeDbl)
			.def("__div__",divNodeMpfr)
			.def("__rdiv__",divNodeMpfr)
			.def("__div__",divNodeInt)
			.def("__rdiv__",divNodeInt)
			.def("__idiv__",&idivNodeNode)
			.def("__idiv__",&idivNodeDouble)
			.def("__idiv__",idivMultNode)
			
			.def("__neg__", negNode)

			.def("__pow__",powNodeNode)
			.def("__pow__",powNodeDouble)
			.def("__pow__",powNodeDbl)
			.def("__pow__",powNodeMpfr)
			.def("__pow__",powNodeInt)
			

			;
			
			def("exp", expNodeNode);
			def("log", logNodeNode);
			def("sin", sinNodeNode);
			def("asin", asinNodeNode);
			def("cos", cosNodeNode);
			def("acos", acosNodeNode);
			def("tan", tanNodeNode);
			def("atan", atanNodeNode);
			
		}
		
		

		
	} // re: python
} // re: bertini
