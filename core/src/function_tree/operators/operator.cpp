//This file is part of Bertini 2.
//
//operator.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//operator.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with operator.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

#include "bertini2/function_tree/operators/operator.hpp"


namespace bertini {
namespace node{	


void UnaryOperator::SetOperand(std::shared_ptr<Node> n)
{
	operand_ = n;
}


//Return the only child for the unary operator
std::shared_ptr<Node> UnaryOperator::Operand() const
{
	return operand_;
}



// The default unary operator is a non-polynomial function of its operand -- sqrt, exp, log,
// the trigonometric functions.  Applied to something free of the variable(s) asked about it
// is a constant (degree 0); applied to anything else it is not a polynomial (degree -1).
// The rule lives here once; the two degree-PRESERVING unaries, NegateOperator and
// IntegerPowerOperator, override.  (The previous group version summed the per-variable
// degrees, which equals the total degree only for a single monomial: NegateOperator inherited
// it and reported -(x^2+y^2) as degree 4, inflating System::DegreeBound() -- the same defect
// #397 found in PowerOperator.)
int UnaryOperator::Degree(std::shared_ptr<Variable> const& v) const
{
	return operand_->Degree(v) == 0 ? 0 : -1;
}

int UnaryOperator::Degree(VariableGroup const& vars) const
{
	return operand_->Degree(vars) == 0 ? 0 : -1;
}

std::vector<int> UnaryOperator::MultiDegreeImpl(VariableGroup const& vars) const
{
	
	std::vector<int> deg(vars.size());
	for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
	{
		*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
	}
	return deg;
}


std::size_t UnaryOperator::HashImpl() const
{
	std::size_t h = typeid(*this).hash_code();   // distinguishes Sin/Cos/Tan/Exp/Log/Sqrt/Negate/...
	HashCombine(h, operand_->Hash());
	return h;
}

bool UnaryOperator::IsSame(Node const& other) const
{
	if (typeid(*this) != typeid(other))
		return false;
	// same concrete unary type -> safe to view as UnaryOperator and compare operand identity
	return operand_.get() == static_cast<UnaryOperator const&>(other).operand_.get();
}


bool UnaryOperator::IsHomogeneous(std::shared_ptr<Variable> const& v) const
{
	if (Degree(v)==0)
	{
		return true;
	}
	else
		return false;
}


/**
Check for homogeneity, with respect to a variable group.
*/
bool UnaryOperator::IsHomogeneous(VariableGroup const& vars) const
{
	if (Degree(vars)==0)
	{
		return true;
	}
	else
		return false;
}



/**
 Change the precision of this variable-precision tree node.
 
 \param prec the number of digits to change precision to.
 */



////////////
//
//  Nary 
//
////////////


// Add an operand onto the container for this operator
void NaryOperator::AddOperand(std::shared_ptr<Node> n)
{
	operands_.push_back(std::move(n));
}






size_t NaryOperator::NumOperands() const
{
	return operands_.size();
}

std::shared_ptr<Node> NaryOperator::FirstOperand() const
{
	return operands_[0];
}




 /**
 Change the precision of this variable-precision tree node.
 
 \param prec the number of digits to change precision to.
 */


void NaryOperator::PrecisionChangeSpecific(unsigned /*prec*/) const
{}


} // re: namespace node	
} // re: bertini namespace
