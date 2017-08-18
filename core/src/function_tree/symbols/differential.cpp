//This file is part of Bertini 2.
//
//src/function_tree/symbols/differential.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/symbols/differential.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/symbols/differential.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame
// Jeb Collins, West Texas A&M


#include "function_tree/symbols/differential.hpp"




namespace bertini{
	namespace node{

Differential::Differential(std::shared_ptr<const Variable> diff_variable, std::string var_name) : NamedSymbol('d'+ var_name), differential_variable_(diff_variable)
{
}

void Differential::Reset() const
{
	Node::ResetStoredValues();
}

const std::shared_ptr<const Variable>& Differential::GetVariable() const 
{
	return differential_variable_;
}	

void Differential::print(std::ostream & target) const
{
	target << name();
}


std::shared_ptr<Node> Differential::Differentiate(std::shared_ptr<Variable> const& v) const
{
	throw std::runtime_error("differentiating a differential is not correctly implemented yet");
	return MakeInteger(0);
}

/**
Compute the degree with respect to a single variable.   For differentials, the degree is 0.
*/
int Differential::Degree(std::shared_ptr<Variable> const& v) const
{
	return 0;
}

int Differential::Degree(VariableGroup const& vars) const
{
	return 0;
}

std::vector<int> Differential::MultiDegree(VariableGroup const& vars) const
{
	return std::vector<int>(vars.size(),0);
}

void Differential::Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar)
{
	
}

bool Differential::IsHomogeneous(std::shared_ptr<Variable> const& v) const
{
	return true;
}

/**
Check for homogeneity, with respect to a variable group.
*/
bool Differential::IsHomogeneous(VariableGroup const& vars) const
{
	return true;
}


/**
 Change the precision of this variable-precision tree node.
 
 \param prec the number of digits to change precision to.
 */
void Differential::precision(unsigned int prec) const
{
	auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
	val_pair.first.precision(prec);
}

// This should never be called for a Differential.  Only for Jacobians.
dbl Differential::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	if(differential_variable_ == diff_variable)
	{
		return 1.0;
	}
	else
	{
		return 0.0;
	}
}

void Differential::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	if(differential_variable_ == diff_variable)
	{
		evaluation_value = 1.0;
	}
	else
	{
		evaluation_value = 0.0;
	}
}


mpfr Differential::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	if(differential_variable_ == diff_variable)
	{
		return mpfr(1);
	}
	else
	{
		return mpfr(0);
	}
}

void Differential::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	if(differential_variable_ == diff_variable)
	{
		evaluation_value.SetOne();
	}
	else
	{
		evaluation_value.SetZero();
	}
}
	} // re: namespace node
} // re: namespace bertini
