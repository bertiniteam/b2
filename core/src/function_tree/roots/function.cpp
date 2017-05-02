//This file is part of Bertini 2.
//
//src/function_tree/roots/jacobian.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/roots/jacobian.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/roots/jacobian.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame
// Jeb Collins, West Texas A&M


#include "function_tree/roots/jacobian.hpp"




namespace bertini{
namespace node{

Function::Function(std::string new_name) : NamedSymbol(new_name)
{ }


Function::Function(const std::shared_ptr<Node> & entry) : NamedSymbol("unnamed_function"), entry_node_(entry)
{
}

const std::shared_ptr<Node> & Function::entry_node() const
{
	return entry_node_;
}


void Function::print(std::ostream & target) const
{
	EnsureNotEmpty();
	entry_node_->print(target);
}


void Function::Reset() const
{
	EnsureNotEmpty();
	
	Node::ResetStoredValues();
	entry_node_->Reset();
}

void Function::SetRoot(std::shared_ptr<Node> const& entry)
{
	entry_node_ = entry;
}

void Function::EnsureNotEmpty() const
{
	if (entry_node_==nullptr)
	{
		throw std::runtime_error("Function node type has empty root node");
	}
}

std::shared_ptr<Node> Function::Differentiate(std::shared_ptr<Variable> const& v) const
{
	return entry_node_->Differentiate(v);
}

/**
Compute the degree of a node.  For functions, the degree is the degree of the entry node.
*/
int Function::Degree(std::shared_ptr<Variable> const& v) const
{
	return entry_node_->Degree(v);
}


int Function::Degree(VariableGroup const& vars) const
{
	return entry_node_->Degree(vars);
}

std::vector<int> Function::MultiDegree(VariableGroup const& vars) const
{
	
	std::vector<int> deg(vars.size());
	for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
	{
		*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
	}
	return deg;
}


void Function::Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar)
{
	entry_node_->Homogenize(vars, homvar);
}

bool Function::IsHomogeneous(std::shared_ptr<Variable> const& v) const
{
	return entry_node_->IsHomogeneous(v);
}

/**
Check for homogeneity, with respect to a variable group.
*/
bool Function::IsHomogeneous(VariableGroup const& vars) const
{
	return entry_node_->IsHomogeneous(vars);
}


/**
 Change the precision of this variable-precision tree node.
 
 \param prec the number of digits to change precision to.
 */
void Function::precision(unsigned int prec) const
{
	auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
	if (val_pair.first.precision()==prec)
		return;
	else{
		val_pair.first.precision(prec);
		entry_node_->precision(prec);
	}
	
}


/**
 Calls FreshEval on the entry node to the tree.
 */
dbl Function::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return entry_node_->Eval<dbl>(diff_variable);
}

/**
 Calls FreshEval in place on the entry node to the tree.
 */
void Function::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	entry_node_->EvalInPlace<dbl>(evaluation_value, diff_variable);
}


/**
 Calls FreshEval on the entry node to the tree.
 */
mpfr Function::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return entry_node_->Eval<mpfr>(diff_variable);
}

/**
 Calls FreshEval in place on the entry node to the tree.
 */
void Function::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	entry_node_->EvalInPlace<mpfr>(evaluation_value, diff_variable);
}



} // re: namespace node
} // re: namespace bertini



