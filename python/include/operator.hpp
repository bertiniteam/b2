// python/function_tree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/function_tree.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/function_tree.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Jeb Collins
//  West Texas A&M University
//  Mathematics
//  Fall 2015
//
//
//  python/operator.hpp:  the header file for the python interface for operator classes.

#ifndef BERTINI_PYTHON_OPERATOR_HPP
#define BERTINI_PYTHON_OPERATOR_HPP


#include "function_tree.hpp"

namespace bertini{
	namespace python{
		
		using namespace boost::python;
		using namespace bertini::node;
		
		// UnaryOperator class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(unOpDeg1Overloads, node::UnaryOperator::Degree, 0, 1)
		int (node::UnaryOperator::*unOpDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::UnaryOperator::Degree;
		int (node::UnaryOperator::*unOpDeg2)(VariableGroup const&) const  = &bertini::node::UnaryOperator::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(unOpIsHom1Overloads, node::UnaryOperator::IsHomogeneous, 0, 1)
		bool (node::UnaryOperator::*unOpIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::UnaryOperator::IsHomogeneous;
		bool (node::UnaryOperator::*unOpIsHom2)(VariableGroup const& vars) const= &bertini::node::UnaryOperator::IsHomogeneous;
		
		// SumOperator class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(sumDeg1Overloads, node::SumOperator::Degree, 0, 1)
		int (node::SumOperator::*sumDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::SumOperator::Degree;
		int (node::SumOperator::*sumDeg2)(VariableGroup const&) const  = &bertini::node::SumOperator::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(sumIsHom1Overloads, node::SumOperator::IsHomogeneous, 0, 1)
		bool (node::SumOperator::*sumIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::SumOperator::IsHomogeneous;
		bool (node::SumOperator::*sumIsHom2)(VariableGroup const& vars) const= &bertini::node::SumOperator::IsHomogeneous;
		void (node::SumOperator::*sumAddChild1)(std::shared_ptr<Node> child) = &node::SumOperator::AddChild;
		void (node::SumOperator::*sumAddChild2)(std::shared_ptr<Node> child, bool) = &node::SumOperator::AddChild;
		
		// NegOperator class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(negIsHom1Overloads, node::NegateOperator::IsHomogeneous, 0, 1)
		bool (node::NegateOperator::*negIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::NegateOperator::IsHomogeneous;
		bool (node::NegateOperator::*negIsHom2)(VariableGroup const& vars) const= &bertini::node::NegateOperator::IsHomogeneous;
		
		// MultOperator class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(multDeg1Overloads, node::MultOperator::Degree, 0, 1)
		int (node::MultOperator::*multDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::MultOperator::Degree;
		int (node::MultOperator::*multDeg2)(VariableGroup const&) const  = &bertini::node::MultOperator::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(multIsHom1Overloads, node::MultOperator::IsHomogeneous, 0, 1)
		bool (node::MultOperator::*multIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::MultOperator::IsHomogeneous;
		bool (node::MultOperator::*multIsHom2)(VariableGroup const& vars) const= &bertini::node::MultOperator::IsHomogeneous;
		void (node::MultOperator::*multAddChild1)(std::shared_ptr<Node> child) = &node::MultOperator::AddChild;
		void (node::MultOperator::*multAddChild2)(std::shared_ptr<Node> child, bool) = &node::MultOperator::AddChild;

		// PowerOperator class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(powDeg1Overloads, node::PowerOperator::Degree, 0, 1)
		int (node::PowerOperator::*powDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::PowerOperator::Degree;
		int (node::PowerOperator::*powDeg2)(VariableGroup const&) const  = &bertini::node::PowerOperator::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(powIsHom1Overloads, node::PowerOperator::IsHomogeneous, 0, 1)
		bool (node::PowerOperator::*powIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::PowerOperator::IsHomogeneous;
		bool (node::PowerOperator::*powIsHom2)(VariableGroup const& vars) const= &bertini::node::PowerOperator::IsHomogeneous;

		// IntegerPowerOperator class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(intpowDeg1Overloads, node::IntegerPowerOperator::Degree, 0, 1)
		int (node::IntegerPowerOperator::*intpowDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::IntegerPowerOperator::Degree;
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(intpowIsHom1Overloads, node::IntegerPowerOperator::IsHomogeneous, 0, 1)
		bool (node::IntegerPowerOperator::*intpowIsHom1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::IntegerPowerOperator::IsHomogeneous;
		bool (node::IntegerPowerOperator::*intpowIsHom2)(VariableGroup const& vars) const= &bertini::node::IntegerPowerOperator::IsHomogeneous;

		// SqrtOperator class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(sqrtDeg1Overloads, node::SqrtOperator::Degree, 0, 1)
		int (node::SqrtOperator::*sqrtDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::SqrtOperator::Degree;

		// ExpOperator class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(expDeg1Overloads, node::ExpOperator::Degree, 0, 1)
		int (node::ExpOperator::*expDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::ExpOperator::Degree;

		// LogOperator class
		BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(logDeg1Overloads, node::LogOperator::Degree, 0, 1)
		int (node::LogOperator::*logDeg1)(std::shared_ptr<node::Variable> const&) const= &bertini::node::LogOperator::Degree;

		
		
		
		
		
		
	} // re: python
} // re: bertini


#endif
