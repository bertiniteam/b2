//This file is part of Bertini 2.
//
//python/node_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/node_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/node_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/node_export.hpp:  Header file for exposing Node class to python.

#ifndef BERTINI_PYTHON_NODE_EXPORT_HPP
#define BERTINI_PYTHON_NODE_EXPORT_HPP

#include <bertini2/function_tree.hpp>
#include <bertini2/function_tree/node.hpp>


#include "python_common.hpp"

namespace bertini{
	namespace python{


		using namespace bertini::node;

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
			static unsigned GetPrecision(NodeBaseT const& self) {return self.precision();}
			static void SetPrecision(NodeBaseT & self, unsigned p){self.precision(p);}

			static Nodeptr Diff0(NodeBaseT& self) { return self.Differentiate();}
			Nodeptr (NodeBaseT::*Diff1)(std::shared_ptr<Variable> const&) const= &NodeBaseT::Differentiate;

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

			// Use templating to return member function pointer to Eval<T>
			template <typename T>
			using Eval1_ptr = T (NodeBaseT::*)(std::shared_ptr<Variable> const&) const;

			template <typename T>
			static Eval1_ptr<T> return_Eval1_ptr()
			{
				return &NodeBaseT::template Eval<T>;
			};



			// Addition operators
			Nodeptr(*addNodeNode)(Nodeptr, const Nodeptr&) = &(operator+);
			Nodeptr(*addNodeMpfr)(Nodeptr, const mpfr&) = &(operator+);
			Nodeptr(*addNodeInt)(Nodeptr, int) = &(operator+);
			static Nodeptr iaddNodeNode(Nodeptr  lhs, const Nodeptr & rhs)
			{
				return lhs += rhs;
			}
			// static Nodeptr iaddNodeDouble(Nodeptr  lhs, double rhs)
			// {
			// 	return lhs += rhs;
			// }
			static SumOperator iaddSumNode(SumOperator  lhs, const Nodeptr & rhs)
			{
				return lhs += rhs;
			}

			// Subtraction operators
			Nodeptr(*subNodeNode)(Nodeptr, const Nodeptr&) = &(operator-);
			Nodeptr(*subNodeMpfr)(Nodeptr, mpfr) = &(operator-);
			Nodeptr(*subNodeInt)(Nodeptr, int) = &(operator-);
			static Nodeptr isubNodeNode(Nodeptr  lhs, const Nodeptr & rhs)
			{
				return lhs -= rhs;
			}

			static SumOperator isubSumNode(SumOperator  lhs, const Nodeptr & rhs)
			{
				return lhs -= rhs;
			}


			// Negate operator
			Nodeptr(*negNode)(const Nodeptr &) = &(operator-);




			// Multiplication operators
			Nodeptr(*multNodeNode)(Nodeptr, const Nodeptr&) = &(operator*);
			Nodeptr(*multNodeMpfr)(Nodeptr, mpfr) = &(operator*);
			Nodeptr(*multNodeInt)(Nodeptr, int) = &(operator*);
			static Nodeptr imultNodeNode(Nodeptr  lhs, const Nodeptr & rhs)
			{
				return lhs *= rhs;
			}

			Nodeptr(*imultMultNode)(std::shared_ptr<node::MultOperator> &, const Nodeptr &) = &(operator*=);

			// Division operators
			Nodeptr(*divNodeNode)(Nodeptr, const Nodeptr&) = &(operator/);
			Nodeptr(*divNodeMpfr)(Nodeptr, mpfr) = &(operator/);
			Nodeptr(*divNodeInt)(Nodeptr, int) = &(operator/);
			static Nodeptr idivNodeNode(Nodeptr  lhs, const Nodeptr & rhs)
			{
				return lhs /= rhs;
			}

			Nodeptr(*idivMultNode)(std::shared_ptr<node::MultOperator> &, const Nodeptr &) = &(operator/=);

			// Power operators
			Nodeptr(*powNodeNode)(const Nodeptr &, const Nodeptr&) = &pow;
			Nodeptr(*powNodeMpfr)(const Nodeptr&, mpfr) = &pow;
			Nodeptr(*powNodeInt)( Nodeptr const&, int) = &pow;

			// Transcendental operators
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
