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
		
		void ExportNode();
		
		
		
		
		
		
		
		
		
		
		template<typename NodeBaseT>
		class NodeVisitor: public def_visitor<NodeVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:
			
			static int Deg0(NodeBaseT& self) { return self.Degree();}
			int (NodeBaseT::*Deg1)(std::shared_ptr<Variable> const&) const= &NodeBaseT::Degree;
			int (NodeBaseT::*Deg2)(VariableGroup const&) const  = &NodeBaseT::Degree;
			
			static bool IsHom0(NodeBaseT& self) { return self.IsHomogeneous();}
			bool (NodeBaseT::*IsHom1)(std::shared_ptr<Variable> const&) const= &NodeBaseT::IsHomogeneous;
			bool (NodeBaseT::*IsHom2)(VariableGroup const& vars) const= &NodeBaseT::IsHomogeneous;
			
			static bool IsPoly0(NodeBaseT& self) { return self.IsPolynomial();}
			bool (NodeBaseT::*IsPoly1)(std::shared_ptr<Variable> const&) const= &NodeBaseT::IsPolynomial;
			bool (NodeBaseT::*IsPoly2)(VariableGroup const& vars) const= &NodeBaseT::IsPolynomial;

			// Can't create member function pointer to Eval with zero arguments because implementation
			// uses default arguments
			template <typename T>
			static T Eval0(NodeBaseT& self) { return self.template Eval<T>();}
			
			template <typename T>
			using Eval1_ptr = T (NodeBaseT::*)(std::shared_ptr<Variable>);
			template <typename T>
			static Eval1_ptr<T> return_Eval1_ptr()
			{
				return &NodeBaseT::template Eval<T>;
			};
			


			// Addition operators
			Nodeptr(*addNodeNode)(Nodeptr, const Nodeptr&) = &operator+;
			Nodeptr(*addNodeDouble)(Nodeptr, double) = &operator+;
			Nodeptr(*addNodeDbl)(Nodeptr, dbl) = &operator+;
			Nodeptr(*addNodeMpfr)(Nodeptr, mpfr) = &operator+;
			Nodeptr(*addNodeInt)(Nodeptr, int) = &operator+;
			static Nodeptr iaddNodeNode(Nodeptr  lhs, const Nodeptr & rhs)
			{
				return lhs += rhs;
			}
			static Nodeptr iaddNodeDouble(Nodeptr  lhs, double rhs)
			{
				return lhs += rhs;
			}
			static SumOperator iaddSumNode(SumOperator  lhs, const Nodeptr & rhs)
			{
				return lhs += rhs;
			}
			
			// Subtraction operators
			Nodeptr(*subNodeNode)(Nodeptr, const Nodeptr&) = &operator-;
			Nodeptr(*subNodeDouble)(Nodeptr, double) = &operator-;
			Nodeptr(*subNodeDbl)(Nodeptr, dbl) = &operator-;
			Nodeptr(*subNodeMpfr)(Nodeptr, mpfr) = &operator-;
			Nodeptr(*subNodeInt)(Nodeptr, int) = &operator-;
			static Nodeptr isubNodeNode(Nodeptr  lhs, const Nodeptr & rhs)
			{
				return lhs -= rhs;
			}
			static Nodeptr isubNodeDouble(Nodeptr  lhs, double rhs)
			{
				return lhs -= rhs;
			}
			static SumOperator isubSumNode(SumOperator  lhs, const Nodeptr & rhs)
			{
				return lhs -= rhs;
			}
			
			
			// Negate operator
			Nodeptr(*negNode)(const Nodeptr &) = &operator-;
			
			
			
			
			// Multiplication operators
			Nodeptr(*multNodeNode)(Nodeptr, const Nodeptr&) = &operator*;
			Nodeptr(*multNodeDouble)(Nodeptr, double) = &operator*;
			Nodeptr(*multNodeDbl)(Nodeptr, dbl) = &operator*;
			Nodeptr(*multNodeMpfr)(Nodeptr, mpfr) = &operator*;
			Nodeptr(*multNodeInt)(Nodeptr, int) = &operator*;
			static Nodeptr imultNodeNode(Nodeptr  lhs, const Nodeptr & rhs)
			{
				return lhs *= rhs;
			}
			static Nodeptr imultNodeDouble(Nodeptr  lhs, double rhs)
			{
				return lhs *= rhs;
			}
			Nodeptr(*imultMultNode)(std::shared_ptr<node::MultOperator> &, const Nodeptr &) = &operator*=;
			
			// Division operators
			Nodeptr(*divNodeNode)(Nodeptr, const Nodeptr&) = &operator/;
			Nodeptr(*divNodeDouble)(Nodeptr, double) = &operator/;
			Nodeptr(*divNodeDbl)(Nodeptr, dbl) = &operator/;
			Nodeptr(*divNodeMpfr)(Nodeptr, mpfr) = &operator/;
			Nodeptr(*divNodeInt)(Nodeptr, int) = &operator/;
			static Nodeptr idivNodeNode(Nodeptr  lhs, const Nodeptr & rhs)
			{
				return lhs /= rhs;
			}
			static Nodeptr idivNodeDouble(Nodeptr  lhs, double rhs)
			{
				return lhs /= rhs;
			}
			Nodeptr(*idivMultNode)(std::shared_ptr<node::MultOperator> &, const Nodeptr &) = &operator/=;
			
			// Power operators
			Nodeptr(*powNodeNode)(const Nodeptr &, const Nodeptr&) = &pow;
			Nodeptr(*powNodeDouble)(const Nodeptr&, double) = &pow;
			Nodeptr(*powNodeDbl)(const Nodeptr&, dbl) = &pow;
			Nodeptr(*powNodeMpfr)(const Nodeptr&, mpfr) = &pow;
			Nodeptr(*powNodeInt)( Nodeptr const&, int) = &pow;
			
			// Transcindental operators
			Nodeptr(*expNodeNode)(const Nodeptr &) = &exp;
			Nodeptr(*logNodeNode)(const Nodeptr &) = &log;
			Nodeptr(*sinNodeNode)(const Nodeptr &) = &sin;
			Nodeptr(*asinNodeNode)(const Nodeptr &) = &asin;
			Nodeptr(*cosNodeNode)(const Nodeptr &) = &cos;
			Nodeptr(*acosNodeNode)(const Nodeptr &) = &acos;
			Nodeptr(*tanNodeNode)(const Nodeptr &) = &tan;
			Nodeptr(*atanNodeNode)(const Nodeptr &) = &atan;
			
			
		};
		
	}
}

		
		
		
		
		
		

#endif
