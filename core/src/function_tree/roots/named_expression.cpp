//This file is part of Bertini 2.
//
//src/function_tree/roots/named_expression.cpp is free software: you can redistribute it and/or
//modify it under the terms of the GNU General Public License as published by the Free Software
//Foundation, either version 3 of the License, or (at your option) any later version.
//
//src/function_tree/roots/named_expression.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
//PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License along with
//src/function_tree/roots/named_expression.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire
// Jeb Collins, West Texas A&M


#include "bertini2/function_tree.hpp"

namespace bertini{
namespace node{


dbl NamedExpression::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return entry_node_->Eval<dbl>(diff_variable);
}

void NamedExpression::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	entry_node_->EvalInPlace<dbl>(evaluation_value, diff_variable);
}


mpfr_complex NamedExpression::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return entry_node_->Eval<mpfr_complex>(diff_variable);
}

void NamedExpression::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	entry_node_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
}


const std::shared_ptr<Node> & NamedExpression::EntryNode() const
{
	assert (entry_node_!=nullptr);

	return entry_node_;
}


void NamedExpression::Reset() const
{
	EnsureNotEmpty();

	Node::ResetStoredValues();
	entry_node_->Reset();
}


void NamedExpression::EnsureNotEmpty() const
{
	if (entry_node_==nullptr)
	{
		throw std::runtime_error("NamedExpression node type has empty entry node");
	}
}

std::shared_ptr<Node> NamedExpression::Differentiate(std::shared_ptr<Variable> const& v) const
{
	return entry_node_->Differentiate(v);
}

int NamedExpression::Degree(std::shared_ptr<Variable> const& v) const
{
	return entry_node_->Degree(v);
}


int NamedExpression::Degree(VariableGroup const& vars) const
{
	return entry_node_->Degree(vars);
}

std::vector<int> NamedExpression::MultiDegree(VariableGroup const& vars) const
{

	std::vector<int> deg(vars.size());
	for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
	{
		*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
	}
	return deg;
}


std::shared_ptr<Node> NamedExpression::Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const
{
	auto homogenized_entry = entry_node_->Homogenized(vars, homvar);
	if (homogenized_entry == entry_node_)
		return std::const_pointer_cast<Node>(shared_from_this());  // unchanged -- share
	// functional: a fresh NamedExpression wrapping the homogenized entry; the original is
	// untouched, so a user (or another System) holding this node never observes it change.
	return NamedExpression::Make(homogenized_entry, name());
}

bool NamedExpression::IsHomogeneous(std::shared_ptr<Variable> const& v) const
{
	return entry_node_->IsHomogeneous(v);
}

bool NamedExpression::IsHomogeneous(VariableGroup const& vars) const
{
	return entry_node_->IsHomogeneous(vars);
}


void NamedExpression::precision(unsigned int prec) const
{
	auto& val_pair = std::get< std::pair<mpfr_complex,bool> >(current_value_);
	if (val_pair.first.precision()==prec)
		return;
	else{
		val_pair.first.precision(prec);
		entry_node_->precision(prec);
	}

}


} // re: namespace node
} // re: namespace bertini
