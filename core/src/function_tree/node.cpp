//This file is part of Bertini 2.
//
//node.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//node.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with node.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame
// Jeb Collins, West Texas A&M


#include "function_tree.hpp"

BOOST_CLASS_EXPORT(bertini::node::Variable)
BOOST_CLASS_EXPORT(bertini::node::Differential)

BOOST_CLASS_EXPORT(bertini::node::Float)
BOOST_CLASS_EXPORT(bertini::node::Integer)
BOOST_CLASS_EXPORT(bertini::node::Rational)

BOOST_CLASS_EXPORT(bertini::node::special_number::Pi)
BOOST_CLASS_EXPORT(bertini::node::special_number::E)

BOOST_CLASS_EXPORT(bertini::node::Function)
BOOST_CLASS_EXPORT(bertini::node::Jacobian)

BOOST_CLASS_EXPORT(bertini::node::SinOperator)
BOOST_CLASS_EXPORT(bertini::node::ArcSinOperator)
BOOST_CLASS_EXPORT(bertini::node::CosOperator)
BOOST_CLASS_EXPORT(bertini::node::ArcCosOperator)
BOOST_CLASS_EXPORT(bertini::node::TanOperator)
BOOST_CLASS_EXPORT(bertini::node::ArcTanOperator)


BOOST_CLASS_EXPORT(bertini::node::SumOperator)
BOOST_CLASS_EXPORT(bertini::node::NegateOperator)
BOOST_CLASS_EXPORT(bertini::node::MultOperator)
BOOST_CLASS_EXPORT(bertini::node::PowerOperator)
BOOST_CLASS_EXPORT(bertini::node::IntegerPowerOperator)
BOOST_CLASS_EXPORT(bertini::node::SqrtOperator)
BOOST_CLASS_EXPORT(bertini::node::ExpOperator)



namespace bertini{
namespace node{

	unsigned Node::ReduceDepth()
	{
		return 0;
	}


	template<typename T>
	void Node::EvalInPlace(T& eval_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		auto& val_pair = std::get< std::pair<T,bool> >(current_value_);
		if(!val_pair.second)
		{
			detail::FreshEvalSelector<T>::RunInPlace(val_pair.first, *this,diff_variable);
			val_pair.second = true;
		}
		eval_value = val_pair.first;
	}

	template void Node::EvalInPlace<dbl>(dbl&, std::shared_ptr<Variable> const&) const;
	template void Node::EvalInPlace<mpfr>(mpfr&, std::shared_ptr<Variable> const&) const;

	unsigned Node::precision() const
	{
		return std::get<std::pair<mpfr,bool> >(current_value_).first.precision();
	}

	bool Node::IsPolynomial(std::shared_ptr<Variable> const&v) const
	{
		return Degree(v)>=0;
	}

	bool Node::IsPolynomial(VariableGroup const&v) const
	{
		return Degree(v)>=0;
	}

	void Node::ResetStoredValues() const
	{
		std::get< std::pair<dbl,bool> >(current_value_).second = false;
		std::get< std::pair<mpfr,bool> >(current_value_).second = false;
	}

	Node::Node()
	{
		std::get<std::pair<dbl,bool> >(current_value_).second = false;
		std::get<std::pair<mpfr,bool> >(current_value_).second = false;
	}
} // namespace node
} // namespace bertini


