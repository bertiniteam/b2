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
// Copyright(C) Bertini2 Development Team
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
//  silviana amethyst
//  UWEC
//  2017, Spring 2018
//
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
			
			std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const { return this->get_override("Homogenized")(vars, homvar); }
			
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
			.def("reset", &NodeBaseT::Reset, (arg("self")),"reset (downward) the values of a function tree so that when the next eval_mp or eval_d is called, the tree re-computes")
			.def("precision", &GetPrecision, (arg("self")),"")
			.def("precision", SetPrecision, (arg("self"), arg("precision")),"")
			.def("degree", &Deg0, (arg("self")),"compute the algebraic degree of node in a function tree, with respect to all variables. returns one integer.  negative is non-algebraic.")
			.def("degree", Deg1, (arg("self"),arg("var")),"compute the algebraic degree of node in a function tree, with respect to a particular variable. returns one integer.  negative is non-algebraic.")
			.def("degree", Deg2, (arg("self"),arg("vars")),"compute the algebraic degree of node in a function tree, with respect to a variable group. returns one integer.  negative is non-algebraic.")
			.def("differentiate", Diff0, (arg("self")),"differentiate a node.  is with respect to all variables.  you get a Jacobian back, which represents derivatives wrt all variables simultaneously.")
			.def("differentiate", Diff1, (arg("self")),"differentiate a node with respect to one variable.  You get a regular old Node in a Function Tree back.")
			.def("multidegree", &NodeBaseT::MultiDegree, (arg("self"),arg("vars")),"Compute an integer vector containing the degrees with respect to the variables in `vars`.  Negative entries indicate non-polynomiality")
			.def("homogenized", &NodeBaseT::Homogenized, (arg("self"),arg("vars"), arg("homvar")), "Return a NEW homogenized copy of this function tree (non-mutating) with respect to the variables in `vars` using the homogenizing variable `homvar`.  Degree-deficient terms are padded with powers of `homvar` so all terms share the same degree.  The original tree is left untouched.")
			.def("is_homogeneous", IsHom0,(arg("self")), "test if this Node is homogeneous with respect to all Variables.")
			.def("is_homogeneous", IsHom1,(arg("self"),arg("var")), "test if this Node is homogeneous with respect to Variable `var`.")
			.def("is_homogeneous", IsHom2,(arg("self"),arg("vars")), "test if this Node is homogeneous with respect to the Variables in `vars`.")
			.def("is_polynomial", IsPoly0,(arg("self")), "test if this Node is polynomial with respect to all Variables.")
			.def("is_polynomial", IsPoly1,(arg("self"),arg("var")), "test if this Node is polynomial with respect to Variable `var`.")
			.def("is_polynomial", IsPoly2,(arg("self"),arg("vars")), "test if this Node is polynomial with respect to Variables `vars`.")

			.def("eval_d", &Eval0<dbl>, (arg("self")), "evaluate in double precision.  uses the values of variables already set in a preceding call to `var.set_current_value()`")
			.def("eval_d", return_Eval1_ptr<dbl>(), (arg("self"),arg("var")), "evaluate the derivative of this node with respect to variable `var` in double precision.  uses the values of variables already set in a preceding call to `var.set_current_value()`")
			.def("eval_mp", &Eval0<mpfr_complex>, (arg("self")), "evaluate in multiple precision.  uses the values of variables already set in a preceding call to `var.set_current_value()`")
			.def("eval_mp", return_Eval1_ptr<mpfr_complex>(), (arg("self"),arg("var")), "evaluate the derivative of this node with respect to variable `var` in multiple precision.  uses the values of variables already set in a preceding call to `var.set_current_value()`")
			
			.def(self_ns::str(self_ns::self))
			.def(self_ns::repr(self_ns::self))
			
			.def("__add__",addNodeNode)
			.def("__add__",addNodeMpfr)
			.def("__radd__",&raddNodeMpfr)
			.def("__add__",addNodeRat)
			.def("__radd__",&raddNodeRat)
			.def("__add__",addNodeInt)
			.def("__radd__",raddNodeInt)
			.def("__iadd__",&NodeVisitor::iaddNodeNode)
			.def("__iadd__", &NodeVisitor::iaddSumNode)
			
			.def("__sub__",subNodeNode)
			.def("__sub__",subNodeMpfr)
			.def("__rsub__",rsubNodeMpfr)
			.def("__sub__",subNodeRat)
			.def("__rsub__",rsubNodeRat)
			.def("__sub__",subNodeInt)
			.def("__rsub__",rsubNodeInt)
			.def("__isub__",&NodeVisitor::isubNodeNode)
			.def("__isub__", &NodeVisitor::isubSumNode)
			
			.def("__mul__",multNodeNode)
			.def("__mul__",multNodeMpfr)
			.def("__rmul__",rmultNodeMpfr)
			.def("__mul__",multNodeRat)
			.def("__rmul__",rmultNodeRat)
			.def("__mul__",multNodeInt)
			.def("__rmul__",rmultNodeInt)
			.def("__imul__",&NodeVisitor::imultNodeNode)
			.def("__imul__",imultMultNode)
			

			.def("__div__",divNodeNode)
			.def("__truediv__",divNodeNode)
			.def("__itruediv__",&NodeVisitor::idivNodeNode)



			
			.def("__div__",divNodeMpfr)
			.def("__truediv__",divNodeMpfr)

			.def("__rdiv__",rdivNodeMpfr)
			.def("__rtruediv__",rdivNodeMpfr)

			.def("__rdiv__",rdivNodeRat)
			.def("__rtruediv__",rdivNodeRat)

			.def("__div__",divNodeInt)
			.def("__truediv__",divNodeInt)

			.def("__rdiv__",rdivNodeInt)
			.def("__rtruediv__",rdivNodeInt)

			.def("__idiv__",&NodeVisitor::idivNodeNode)
			.def("__itruediv__",&NodeVisitor::idivNodeNode)
			
			.def("__idiv__",idivMultNode)
			.def("__itruediv__",idivMultNode)
			
			.def("__neg__", negNode)
			
			.def("__pow__",powNodeNode)
			.def("__pow__",powNodeMpfr)
			.def("__pow__",powNodeRat)
			.def("__pow__",powNodeInt)
			;
			
			
			def("exp", expNodeNode, "the symbolic exponential operator");
			def("log", logNodeNode, "the symbolic natural log operator");
			def("sin", sinNodeNode, "the symbolic sine operator");
			def("asin", asinNodeNode, "the symbolic arcsine operator");
			def("cos", cosNodeNode, "the symbolic cosine operator");
			def("acos", acosNodeNode, "the symbolic arccosine operator");
			def("tan", tanNodeNode, "the symbolic tangent operator");
			def("atan", atanNodeNode, "the symbolic arctangent operator");
			
		}

		
		
		
		void ExportNode()
		{			
			class_<NodeWrap, boost::noncopyable, Nodeptr >("AbstractNode", no_init)
			.def(NodeVisitor<Node>())
			;
		};
		
		
	} //namespace python
} // namespace bertini
