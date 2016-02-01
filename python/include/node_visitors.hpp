//
//  node_visitors.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 1/30/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_node_visitors_hpp
#define Xcode_b2_node_visitors_hpp

#include <bertini2/function_tree.hpp>
#include <bertini2/function_tree/node.hpp>


#include "python_common.hpp"

namespace bertini{
	namespace python{
		

		using namespace bertini::node;
		
		using Node = Node;
		using Nodeptr = std::shared_ptr<Node>;
		
		
		
		
		struct NodeWrap : Node, wrapper<Node>
		{
			void Reset()
			{
				if (override Reset = this->get_override("Reset"))
					Reset(); // *note*

				Node::Reset();
			}
			
//			void Reset()
//			{
//				if (override Reset = this->get_override("Reset"))
//					Reset(); // *note*
//				
//				Node::Reset();
//			}

			void default_Reset()
			{
				return this->Node::Reset();
			}

			void print(std::ostream& target) const
			{
				this->get_override("print")();
			}
		}; // re: NodeWrap

		
		
		
		
		template<typename NodeBaseT>
		class NodeBaseVisitor: public def_visitor<NodeBaseVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
			{
				cl
				.def("degree", &NodeBaseVisitor::Deg1)
				.def("degree", &NodeBaseVisitor::Deg2)
				.def("reset", &NodeBaseT::Reset)
				.def("differentiate", &NodeBaseT::Differentiate)
				.def("multidegree", &NodeBaseT::MultiDegree)
				.def("homogenize", &NodeBaseT::Homogenize)
				.def("is_homogeneous", &NodeBaseVisitor::IsHom1)
				.def("is_homogeneous", &NodeBaseVisitor::IsHom2)
				.def("precision", &NodeBaseT::precision)
				.def("is_polynomial", &NodeBaseVisitor::IsPoly1)
				.def("is_polynomial", &NodeBaseVisitor::IsPoly2)
				;
			}
			
			
			
			
		private:
			bool Deg1(const object& obj, std::shared_ptr<Variable> const&v = nullptr) const
			{
				const NodeBaseT& self = extract<NodeBaseT>(obj)();
				return self.Degree(v);
			}
			bool Deg2(const object& obj, VariableGroup const& v) const
			{
				const NodeBaseT& self = extract<NodeBaseT>(obj)();
				return self.Degree(v);
			}

			bool IsHom1(const object& obj, std::shared_ptr<Variable> const&v = nullptr) const
			{
				const NodeBaseT& self = extract<NodeBaseT>(obj)();
				return self.IsHomogeneous(v);
			}
			bool IsHom2(const object& obj, VariableGroup const& v) const
			{
				const NodeBaseT& self = extract<NodeBaseT>(obj)();
				return self.IsHomogeneous(v);
			}

			bool IsPoly1(const object& obj, std::shared_ptr<Variable> const&v = nullptr) const
			{
				const NodeBaseT& self = extract<NodeBaseT>(obj)();
				return self.IsPolynomial(v);
			}
			bool IsPoly2(const object& obj, VariableGroup const& v) const
			{
				const NodeBaseT& self = extract<NodeBaseT>(obj)();
				return self.IsPolynomial(v);
			}

//			BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Deg1Overloads, NodeBaseT::Degree, 0, 1)
//			int (NodeBaseT::*Deg1)(std::shared_ptr<Variable> const&) const= &NodeBaseT::Degree;
//			int (NodeBaseT::*Deg2)(VariableGroup const&) const  = &NodeBaseT::Degree;
//			
//			BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(IsHom1Overloads, NodeBaseT::IsHomogeneous, 0, 1)
//			bool (NodeBaseT::*IsHom1)(std::shared_ptr<Variable> const&) const= &NodeBaseT::IsHomogeneous;
//			bool (NodeBaseT::*IsHom2)(VariableGroup const& vars) const= &NodeBaseT::IsHomogeneous;

//			BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(IsPoly1Overloads, NodeBaseT::IsPolynomial, 0, 1)
//			bool (NodeBaseT::*IsPoly1)(std::shared_ptr<Variable> const&) const= &NodeBaseT::IsPolynomial;
//			bool (Node::*IsPoly2)(VariableGroup const& vars) const= &NodeBaseT::IsPolynomial;

		};
		
		
		
		
		
		
		
		template<typename NodeBaseT>
		class NodeVisitor: public def_visitor<NodeVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const
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
			
		private:
			BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NodeEvalOverloadsDbl, Node::template Eval<dbl>, 0, 1);
			BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NodeEvalOverloadsMPFR, Node::template Eval<mpfr>, 0, 1);

			// Addition operators
			std::shared_ptr<Node>(*addNodeNode)(std::shared_ptr<Node>, const std::shared_ptr<Node>&) = &operator+;
			std::shared_ptr<Node>(*addNodeDouble)(std::shared_ptr<Node>, double) = &operator+;
			std::shared_ptr<Node>(*addNodeDbl)(std::shared_ptr<Node>, dbl) = &operator+;
			std::shared_ptr<Node>(*addNodeMpfr)(std::shared_ptr<Node>, mpfr) = &operator+;
			std::shared_ptr<Node>(*addNodeInt)(std::shared_ptr<Node>, int) = &operator+;
			inline std::shared_ptr<Node> iaddNodeNode(std::shared_ptr<Node>  lhs, const std::shared_ptr<Node> & rhs)
			{
				return lhs += rhs;
			}
			inline std::shared_ptr<Node> iaddNodeDouble(std::shared_ptr<Node>  lhs, double rhs)
			{
				return lhs += rhs;
			}
			inline SumOperator iaddSumNode(SumOperator  lhs, const std::shared_ptr<Node> & rhs)
			{
				return lhs += rhs;
			}
			
			// Subtraction operators
			std::shared_ptr<Node>(*subNodeNode)(std::shared_ptr<Node>, const std::shared_ptr<Node>&) = &operator-;
			std::shared_ptr<Node>(*subNodeDouble)(std::shared_ptr<Node>, double) = &operator-;
			std::shared_ptr<Node>(*subNodeDbl)(std::shared_ptr<Node>, dbl) = &operator-;
			std::shared_ptr<Node>(*subNodeMpfr)(std::shared_ptr<Node>, mpfr) = &operator-;
			std::shared_ptr<Node>(*subNodeInt)(std::shared_ptr<Node>, int) = &operator-;
			inline std::shared_ptr<Node> isubNodeNode(std::shared_ptr<Node>  lhs, const std::shared_ptr<Node> & rhs)
			{
				return lhs -= rhs;
			}
			inline std::shared_ptr<Node> isubNodeDouble(std::shared_ptr<Node>  lhs, double rhs)
			{
				return lhs -= rhs;
			}
			inline SumOperator isubSumNode(SumOperator  lhs, const std::shared_ptr<Node> & rhs)
			{
				return lhs -= rhs;
			}
			
			
			// Negate operator
			std::shared_ptr<Node>(*negNode)(const std::shared_ptr<Node> &) = &operator-;
			
			
			
			
			// Multiplication operators
			std::shared_ptr<Node>(*multNodeNode)(std::shared_ptr<Node>, const std::shared_ptr<Node>&) = &operator*;
			std::shared_ptr<Node>(*multNodeDouble)(std::shared_ptr<Node>, double) = &operator*;
			std::shared_ptr<Node>(*multNodeDbl)(std::shared_ptr<Node>, dbl) = &operator*;
			std::shared_ptr<Node>(*multNodeMpfr)(std::shared_ptr<Node>, mpfr) = &operator*;
			std::shared_ptr<Node>(*multNodeInt)(std::shared_ptr<Node>, int) = &operator*;
			inline std::shared_ptr<Node> imultNodeNode(std::shared_ptr<Node>  lhs, const std::shared_ptr<Node> & rhs)
			{
				return lhs *= rhs;
			}
			inline std::shared_ptr<Node> imultNodeDouble(std::shared_ptr<Node>  lhs, double rhs)
			{
				return lhs *= rhs;
			}
			std::shared_ptr<Node>(*imultMultNode)(std::shared_ptr<node::MultOperator> &, const std::shared_ptr<Node> &) = &operator*=;
			
			// Division operators
			std::shared_ptr<Node>(*divNodeNode)(std::shared_ptr<Node>, const std::shared_ptr<Node>&) = &operator/;
			std::shared_ptr<Node>(*divNodeDouble)(std::shared_ptr<Node>, double) = &operator/;
			std::shared_ptr<Node>(*divNodeDbl)(std::shared_ptr<Node>, dbl) = &operator/;
			std::shared_ptr<Node>(*divNodeMpfr)(std::shared_ptr<Node>, mpfr) = &operator/;
			std::shared_ptr<Node>(*divNodeInt)(std::shared_ptr<Node>, int) = &operator/;
			inline std::shared_ptr<Node> idivNodeNode(std::shared_ptr<Node>  lhs, const std::shared_ptr<Node> & rhs)
			{
				return lhs /= rhs;
			}
			inline std::shared_ptr<Node> idivNodeDouble(std::shared_ptr<Node>  lhs, double rhs)
			{
				return lhs /= rhs;
			}
			std::shared_ptr<Node>(*idivMultNode)(std::shared_ptr<node::MultOperator> &, const std::shared_ptr<Node> &) = &operator/=;
			
			// Power operators
			std::shared_ptr<Node>(*powNodeNode)(const std::shared_ptr<Node> &, const std::shared_ptr<Node>&) = &pow;
			std::shared_ptr<Node>(*powNodeDouble)(const std::shared_ptr<Node>&, double) = &pow;
			std::shared_ptr<Node>(*powNodeDbl)(const std::shared_ptr<Node>&, dbl) = &pow;
			std::shared_ptr<Node>(*powNodeMpfr)(const std::shared_ptr<Node>&, mpfr) = &pow;
			std::shared_ptr<Node>(*powNodeInt)( std::shared_ptr<Node> const&, int) = &pow;
			
			// Transcindental operators
			std::shared_ptr<Node>(*expNodeNode)(const std::shared_ptr<Node> &) = &exp;
			std::shared_ptr<Node>(*logNodeNode)(const std::shared_ptr<Node> &) = &log;
			std::shared_ptr<Node>(*sinNodeNode)(const std::shared_ptr<Node> &) = &sin;
			std::shared_ptr<Node>(*asinNodeNode)(const std::shared_ptr<Node> &) = &asin;
			std::shared_ptr<Node>(*cosNodeNode)(const std::shared_ptr<Node> &) = &cos;
			std::shared_ptr<Node>(*acosNodeNode)(const std::shared_ptr<Node> &) = &acos;
			std::shared_ptr<Node>(*tanNodeNode)(const std::shared_ptr<Node> &) = &tan;
			std::shared_ptr<Node>(*atanNodeNode)(const std::shared_ptr<Node> &) = &atan;
			
			
		};
		
	}
}

		
		
		
		
		
		

#endif
