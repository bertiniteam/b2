//
//  node_export.cpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/2/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#include <stdio.h>


#include "node_visitors.hpp"

//#include <bertini2/function_tree/node.hpp>


namespace bertini{
	namespace python{
		
		struct NodeWrap : Node, wrapper<Node>
		{
//			void Reset()
//			{
//				if (override Reset = this->get_override("Reset"))
//					Reset(); // *note*
//				
//				Node::Reset();
//			}
			
			void precision(unsigned int prec)
			{
				this->get_override("precision")(prec);
			}
			
//			void default_Reset()
//			{
//				return this->Node::Reset();
//			}
//			
//			void print(std::ostream& target) const
//			{
//				this->get_override("print")();
//			}
		}; // re: NodeWrap





		template<typename NodeBaseT>
		template<class PyClass>
		void NodeBaseVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			
			//				.def("degree", Deg1, Deg1Overloads())
			//				.def("degree", Deg2)
			.def("reset", &NodeBaseT::Reset)
			//				.def("differentiate", &NodeBaseT::Differentiate)
			//				.def("multidegree", &NodeBaseT::MultiDegree)
			//				.def("homogenize", &NodeBaseT::Homogenize)
			//				.def("is_homogeneous", &NodeBaseVisitor::IsHom1)
			//				.def("is_homogeneous", &NodeBaseVisitor::IsHom2)
			//				.def("is_polynomial", &NodeBaseVisitor::IsPoly1)
			//				.def("is_polynomial", &NodeBaseVisitor::IsPoly2)
			;
		}
		
	
		
		
		
		template<typename NodeBaseT>
		template<class PyClass>
		void NodeVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			NodeBaseVisitor<NodeBaseT>().visit(cl);
			
			cl
			.def("evald", &Node::Eval<dbl>,
				 NodeEvalOverloadsDbl(args("DiffVar"), "Eval's docstring"))
			.def("evalmp", &Node::Eval<mpfr>,
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
			.def("__iadd__",&NodeVisitor::iaddNodeNode)
			.def("__iadd__",&NodeVisitor::iaddNodeDouble)
			.def("__iadd__", &NodeVisitor::iaddSumNode)
			
			.def("__sub__",subNodeNode)
			.def("__sub__",subNodeDouble)
			.def("__rsub__",subNodeDouble)
			.def("__sub__",subNodeDbl)
			.def("__rsub__",subNodeDbl)
			.def("__sub__",subNodeMpfr)
			.def("__rsub__",subNodeMpfr)
			.def("__sub__",subNodeInt)
			.def("__rsub__",subNodeInt)
			.def("__isub__",&NodeVisitor::isubNodeNode)
			.def("__isub__",&NodeVisitor::isubNodeDouble)
			.def("__isub__", &NodeVisitor::isubSumNode)
			
			.def("__mul__",multNodeNode)
			.def("__mul__",multNodeDouble)
			.def("__rmul__",multNodeDouble)
			.def("__mul__",multNodeDbl)
			.def("__rmul__",multNodeDbl)
			.def("__mul__",multNodeMpfr)
			.def("__rmul__",multNodeMpfr)
			.def("__mul__",multNodeInt)
			.def("__rmul__",multNodeInt)
			.def("__imul__",&NodeVisitor::imultNodeNode)
			.def("__imul__",&NodeVisitor::imultNodeDouble)
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
			.def("__idiv__",&NodeVisitor::idivNodeNode)
			.def("__idiv__",&NodeVisitor::idivNodeDouble)
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

		
		
		
		void ExportNode()
		{
			class_<NodeWrap, boost::noncopyable, Nodeptr >("Node", no_init)
//			.def(NodeVisitor<NodeWrap>())
			.def("precision", pure_virtual(&Node::precision) )
			;
			
		};
		
		
	} //namespace python
} // namespace bertini
