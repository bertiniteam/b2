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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of notre dame
// Jeb Collins, West Texas A&M


#include "bertini2/function_tree/symbols/variable.hpp"

#include "bertini2/eigen_extensions.hpp"
#include "bertini2/naming.hpp"



namespace bertini{
	namespace node{
		using ::pow;

// The single funnel for a named Variable (Make("...") reaches here): reject a name
// that is not a well-formed identifier, so an expression/operator/whitespace/leading-
// digit string can never become a variable.  The default ctor's placeholder and
// serialization (which restores the name directly) do not pass through here.
Variable::Variable(std::string new_name) : NamedSymbol(new_name)
{
	ThrowIfInvalidVariableName(new_name);
}

Variable::Variable() : NamedSymbol("unnamed_variable_be_scared")
{ }

std::shared_ptr<Node> Variable::Differentiate(std::shared_ptr<Variable> const& v) const
{
	if (v==nullptr)
		return Differential::Make(std::static_pointer_cast<Variable const>(shared_from_this()), name());
	else
		return v.get() == this ? Integer::Make(1) : Integer::Make(0);
}

std::shared_ptr<Node> Variable::Subs(SubstitutionMap const& substitutions) const
{
	auto it = substitutions.find(name());
	if (it != substitutions.end())
		return it->second;
	return std::const_pointer_cast<Node>(shared_from_this());
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
	for (const auto& iter : vars)
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



bool Variable::IsHomogeneous(std::shared_ptr<Variable> const& /*v*/) const
{
	return true;
}

/**
Check for homogeneity, with respect to a variable group.
*/
bool Variable::IsHomogeneous(VariableGroup const& /*vars*/) const
{
	return true;
}


	} // re: namespace node
} // re: namespace bertini
