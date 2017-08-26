//This file is part of Bertini 2.
//
//python/node_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/node_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/node_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  Dani Brake
//  University of Notre Dame
//
//  python/node_export.cpp:  Source file for exposing Node class to python.

#include <stdio.h>


#include "node_export.hpp"


namespace bertini{
	namespace python{
		
		// Wrapper struct to allow derived classes to overide methods in python
		struct NodeWrap : Node, wrapper<Node>
		{
			void Reset()
			{
				if (override Reset = this->get_override("Reset"))
					Reset(); 
				
				Node::Reset();
			}
			void default_Reset(){ return this->Node::Reset();}
			
			void precision(unsigned int prec) { this->get_override("precision")(prec); }
			
			int Degree(std::shared_ptr<Variable> const& v = nullptr) const {return this->get_override("Degree")(v); }
			int Degree(VariableGroup const& vars) const {return this->get_override("Degree")(vars); }
			
			std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const {return this->get_override("Differentiate")(v); }
			
			std::vector<int> MultiDegree(VariableGroup const& vars) const {return this->get_override("MultiDegree")(vars); }
			
			void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) { this->get_override("Homogenize")(vars, homvar); }
			
			bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const {return this->get_override("IsHomogeneous")(v); }
			bool IsHomogeneous(VariableGroup const& vars) const {return this->get_override("IsHomogeneous")(vars); }
			
			bool IsPolynomial(std::shared_ptr<Variable> const&v = nullptr) const {return this->get_override("IsPolynomial")(v); }
			bool IsPolynomial(VariableGroup const&v) const {return this->get_override("IsPolynomial")(v); }
			

			
		}; // re: NodeWrap





		
	
		
		
		
		template<typename NodeBaseT>
		template<class PyClass>
		void NodeVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("reset", &NodeBaseT::Reset)
			.def("precision", &GetPrecision )
			.def("precision", SetPrecision )
			.def("degree", &Deg0 )
			.def("degree", Deg1)
			.def("degree", Deg2 )
			.def("differentiate", Diff0 )
			.def("differentiate", Diff1 )
			.def("multidegree", &NodeBaseT::MultiDegree )
			.def("homogenize", &NodeBaseT::Homogenize )
			.def("is_homogeneous", IsHom0 )
			.def("is_homogeneous", IsHom1 )
			.def("is_homogeneous", IsHom2 )
			.def("is_polynomial", IsPoly0 )
			.def("is_polynomial", IsPoly1 )
			.def("is_polynomial", IsPoly2 )

			.def("evald", &Eval0<dbl> )
			.def("evald", return_Eval1_ptr<dbl>() )
			.def("evalmp", &Eval0<mpfr> )
			.def("evalmp", return_Eval1_ptr<mpfr>() )
			
			.def(self_ns::str(self_ns::self))
			.def(self_ns::repr(self_ns::self))
			
			.def("__add__",addNodeNode)
			.def("__add__",addNodeMpfr)
			.def("__radd__",addNodeMpfr)
			.def("__add__",addNodeInt)
			.def("__radd__",addNodeInt)
			.def("__iadd__",&NodeVisitor::iaddNodeNode)
			.def("__iadd__", &NodeVisitor::iaddSumNode)
			
			.def("__sub__",subNodeNode)
			.def("__sub__",subNodeMpfr)
			.def("__rsub__",subNodeMpfr)
			.def("__sub__",subNodeInt)
			.def("__rsub__",subNodeInt)
			.def("__isub__",&NodeVisitor::isubNodeNode)
			.def("__isub__", &NodeVisitor::isubSumNode)
			
			.def("__mul__",multNodeNode)
			.def("__mul__",multNodeMpfr)
			.def("__rmul__",multNodeMpfr)
			.def("__mul__",multNodeInt)
			.def("__rmul__",multNodeInt)
			.def("__imul__",&NodeVisitor::imultNodeNode)
			.def("__imul__",imultMultNode)
			
			.def("__div__",divNodeNode)
			.def("__div__",divNodeMpfr)
			.def("__rdiv__",divNodeMpfr)
			.def("__div__",divNodeInt)
			.def("__rdiv__",divNodeInt)
			.def("__idiv__",&NodeVisitor::idivNodeNode)
			.def("__idiv__",idivMultNode)
			
			.def("__neg__", negNode)
			
			.def("__pow__",powNodeNode)
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

		
		
		
		void ExportNode()
		{			
			class_<NodeWrap, boost::noncopyable, Nodeptr >("AbstractNode", no_init)
			.def(NodeVisitor<Node>())
			;
		};
		
		
	} //namespace python
} // namespace bertini
