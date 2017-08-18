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

Jacobian::Jacobian(const std::shared_ptr<Node> & entry) : Function(entry)
{
}


// Evaluate the node.  If flag false, just return value, if flag true
//  run the specific FreshEval of the node, then set flag to false.
template<typename T>
T Jacobian::EvalJ(std::shared_ptr<Variable> const& diff_variable) const
{
		auto& val_pair = std::get< std::pair<T,bool> >(current_value_);

		if(diff_variable == current_diff_variable_ && val_pair.second)
			return val_pair.first;
		else
		{
			current_diff_variable_ = diff_variable;
			Reset();
			detail::FreshEvalSelector<T>::RunInPlace(val_pair.first, *this, diff_variable);
			val_pair.second = true;
			return val_pair.first;
		}						
}

template dbl Jacobian::EvalJ(std::shared_ptr<Variable> const& diff_variable) const;
template mpfr Jacobian::EvalJ(std::shared_ptr<Variable> const& diff_variable) const;



template<typename T>
void Jacobian::EvalJInPlace(T& eval_value, std::shared_ptr<Variable> const& diff_variable) const
{
		auto& val_pair = std::get< std::pair<T,bool> >(current_value_);

		if(diff_variable == current_diff_variable_ && val_pair.second)
			eval_value = val_pair.first;
		else
		{
			current_diff_variable_ = diff_variable;
			Reset();
			detail::FreshEvalSelector<T>::RunInPlace(val_pair.first,*this,diff_variable);
			val_pair.second = true;
			eval_value = val_pair.first;
		}						
}

template void Jacobian::EvalJInPlace( dbl&, std::shared_ptr<Variable> const& diff_variable) const;
template void Jacobian::EvalJInPlace( mpfr&, std::shared_ptr<Variable> const& diff_variable) const;


void Jacobian::Reset() const
{
	EnsureNotEmpty();
	
	Node::ResetStoredValues();
	entry_node_->Reset();
}



	} // re: namespace node
} // re: namespace bertini



