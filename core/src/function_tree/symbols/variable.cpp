//This file is part of Bertini 2.
//
//src/function_tree/symbols/variable.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/symbols/variable.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/symbols/variable.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame
// Jeb Collins, West Texas A&M


#include "function_tree/symbols/variable.hpp"




namespace bertini{
	namespace node{
		using ::pow;
		
Variable::Variable(std::string new_name) : NamedSymbol(new_name)
{ }

template <typename T>
void Variable::set_current_value(T const& val)
{
	assert(Precision(std::get< std::pair<T,bool> >(current_value_).first)==Precision(val) && "precision of value setting into variable doesn't match precision of variable.  is default precision correct?");
	
	std::get< std::pair<T,bool> >(current_value_).first = val;
	std::get< std::pair<T,bool> >(current_value_).second = false;
}

template void Variable::set_current_value<dbl>(dbl const&);
template void Variable::set_current_value<mpfr>(mpfr const&);

std::shared_ptr<Node> Variable::Differentiate() const
{
	return MakeDifferential(shared_from_this(), name());
}

std::shared_ptr<Node> Variable::Derivative(std::shared_ptr<Variable> const& v) const
{
	return v.get() == this ? MakeInteger(1) : MakeInteger(0);
}

void Variable::Reset() const
{
	Node::ResetStoredValues();
}



int Variable::Degree(std::shared_ptr<Variable> const& v) const
{
	if (v)
	{
		if (this == v.get())
			return 1;
		else
			return 0;
	}
	else
		return 1;
	
}


int Variable::Degree(VariableGroup const& vars) const
{
	for (auto iter : vars)
		if (this==iter.get())
			return 1;
		
	return 0;
}


std::vector<int> Variable::MultiDegree(VariableGroup const& vars) const
{
	std::vector<int> deg;
	for (auto iter=vars.begin(); iter!=vars.end(); iter++)
		if (this==(*iter).get())
			deg.push_back(1);
		else
			deg.push_back(0);
	return deg;
}


void Variable::Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar)
{
	
}

bool Variable::IsHomogeneous(std::shared_ptr<Variable> const& v) const
{
	return true;
}

/**
Check for homogeneity, with respect to a variable group.
*/
bool Variable::IsHomogeneous(VariableGroup const& vars) const
{
	return true;
}


/**
 Change the precision of this variable-precision tree node.
 
 \param prec the number of digits to change precision to.
 */
void Variable::precision(unsigned int prec) const
{
	auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
	val_pair.first.precision(prec);
}


// Return current value of the variable.
dbl Variable::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return std::get< std::pair<dbl,bool> >(current_value_).first;
}

void Variable::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = std::get< std::pair<dbl,bool> >(current_value_).first;
}


mpfr Variable::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return std::get< std::pair<mpfr,bool> >(current_value_).first;
}

void Variable::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = std::get< std::pair<mpfr,bool> >(current_value_).first;
}



	} // re: namespace node
} // re: namespace bertini
