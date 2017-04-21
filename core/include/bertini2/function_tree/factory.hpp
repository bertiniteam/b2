//This file is part of Bertini 2.
//
//include/bertini2/function_tree/factory.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/function_tree/factory.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/function_tree/factory.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Dani Brake
// University of Notre Dame


/**
\file include/bertini2/function_tree/factory.hpp

\brief Methods for producing nodes of various types. 

These wrapper functions are intended to be re-written to use something like a workspace, sometime in the future.  This will allow one to refer to things like the `1` node, or the `pi` node, or the `f2` node.  Other desired functionality is things like inspecting all declared symbols or variables, or having multiple workspaces.

*/

#pragma once

namespace bertini {

	namespace node{ 
		class Variable;
		class Integer;
		class Float;
		class Rational;
		class Function;
		class Jacobian;
		class Differential;
	}

	/**
	\brief Make a variable.

	\return A node for you!

	\param name The name of the variable to make.  In principle can be anything, but stick to Bertini naming rules.  Don't be silly, you'll cause serialization-deserialization round tripping problems.

	\return A shared pointer to a variable, the name of which you stated.
	*/
	template<typename T>
	std::shared_ptr<node::Variable> MakeVariable(T const& t)
	{
		return std::make_shared<node::Variable>(t);
	}


	/**
	\brief Make an integer.

	\return A node for you!

	\param name The integer to make.
	*/
	template<typename T>
	std::shared_ptr<node::Integer> MakeInteger(T const& t)
	{
		return std::make_shared<node::Integer>(t);
	}


	template<typename ... T>
	std::shared_ptr<node::Rational> MakeRational(T const& ...t)
	{
		return std::make_shared<node::Rational>(t...);
	}


	template<typename ... T>
	std::shared_ptr<node::Float> MakeFloat(T const& ...t)
	{
		return std::make_shared<node::Float>(t...);
	}

	template<typename ... T>
	std::shared_ptr<node::Function> MakeFunction(T const& ...t)
	{
		return std::make_shared<node::Function>(t...);
	}

	template<typename ... T>
	std::shared_ptr<node::Differential> MakeDifferential(T const& ...t)
	{
		return std::make_shared<node::Differential>(t...);
	}

	template<typename ... T>
	std::shared_ptr<node::Jacobian> MakeJacobian(T const& ...t)
	{
		return std::make_shared<node::Jacobian>(t...);
	}

}// namespace bertini
