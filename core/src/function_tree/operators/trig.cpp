//This file is part of Bertini 2.
//
//trig.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//trig.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with trig.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire
// Jeb Collins, West Texas A&M


#include "bertini2/function_tree/operators/trig.hpp"


namespace bertini {
namespace node{	

	int TrigOperator::Degree(std::shared_ptr<Variable> const& v) const
	{
		if (operand_->Degree(v)==0)
		{
			return 0;
		}
		else
		{
			return -1;
		}
	}


	void SinOperator::print(std::ostream & target) const
	{
		target << "sin(";
		operand_->print(target);
		target << ")";
	}

	
	std::shared_ptr<Node> SinOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return SimplifiedMult({
			{cos(operand_), true},
			{operand_->Differentiate(v), true}
		});
	}

	
	
	
	

	void ArcSinOperator::print(std::ostream & target) const
	{
		target << "asin(";
		operand_->print(target);
		target << ")";
	}

	std::shared_ptr<Node> ArcSinOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return SimplifiedMult({
			{operand_->Differentiate(v), true},
			{sqrt(1-pow(operand_,2)), false}
		});
	}


	

	
	

	void CosOperator::print(std::ostream & target) const
	{
		target << "cos(";
		operand_->print(target);
		target << ")";
	}
	
	std::shared_ptr<Node> CosOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return SimplifiedNegate(SimplifiedMult({
			{sin(operand_), true},
			{operand_->Differentiate(v), true}
		}));
	}

	
	
	
	


	void ArcCosOperator::print(std::ostream & target) const
	{
		target << "acos(";
		operand_->print(target);
		target << ")";
	}


	std::shared_ptr<Node> ArcCosOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return SimplifiedNegate(SimplifiedMult({
			{operand_->Differentiate(v), true},
			{sqrt(1-pow(operand_,2)), false}
		}));
	}

	
	
	
	


	void TanOperator::print(std::ostream & target) const
	{
		target << "tan(";
		operand_->print(target);
		target << ")";
	}

	std::shared_ptr<Node> TanOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return SimplifiedMult({
			{operand_->Differentiate(v), true},
			{pow(cos(operand_),2), false}
		});
	}

	
	
	
	


	void ArcTanOperator::print(std::ostream & target) const
	{
		target << "atan(";
		operand_->print(target);
		target << ")";
	}
		

	std::shared_ptr<Node> ArcTanOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return SimplifiedMult({
			{operand_->Differentiate(v), true},
			{1 + pow(operand_,2), false}
		});
	}

	
	
	
	

	// ---- functional Simplified(): rebuild same op type from the simplified operand ----
	std::shared_ptr<Node> SinOperator::Simplified()    const { return SinOperator::Make(operand_->Simplified()); }
	std::shared_ptr<Node> ArcSinOperator::Simplified() const { return ArcSinOperator::Make(operand_->Simplified()); }
	std::shared_ptr<Node> CosOperator::Simplified()    const { return CosOperator::Make(operand_->Simplified()); }
	std::shared_ptr<Node> ArcCosOperator::Simplified() const { return ArcCosOperator::Make(operand_->Simplified()); }
	std::shared_ptr<Node> TanOperator::Simplified()    const { return TanOperator::Make(operand_->Simplified()); }
	std::shared_ptr<Node> ArcTanOperator::Simplified() const { return ArcTanOperator::Make(operand_->Simplified()); }

	// ---- functional Subs(): rebuild same op type from the substituted operand ----
	std::shared_ptr<Node> SinOperator::Subs(SubstitutionMap const& m)    const { return SinOperator::Make(operand_->Subs(m)); }
	std::shared_ptr<Node> ArcSinOperator::Subs(SubstitutionMap const& m) const { return ArcSinOperator::Make(operand_->Subs(m)); }
	std::shared_ptr<Node> CosOperator::Subs(SubstitutionMap const& m)    const { return CosOperator::Make(operand_->Subs(m)); }
	std::shared_ptr<Node> ArcCosOperator::Subs(SubstitutionMap const& m) const { return ArcCosOperator::Make(operand_->Subs(m)); }
	std::shared_ptr<Node> TanOperator::Subs(SubstitutionMap const& m)    const { return TanOperator::Make(operand_->Subs(m)); }
	std::shared_ptr<Node> ArcTanOperator::Subs(SubstitutionMap const& m) const { return ArcTanOperator::Make(operand_->Subs(m)); }

	// ---- functional Homogenized(): rebuild same op type from the homogenized operand ----
	std::shared_ptr<Node> SinOperator::Homogenized(VariableGroup const& v, std::shared_ptr<Variable> const& h)    const { return SinOperator::Make(operand_->Homogenized(v,h)); }
	std::shared_ptr<Node> ArcSinOperator::Homogenized(VariableGroup const& v, std::shared_ptr<Variable> const& h) const { return ArcSinOperator::Make(operand_->Homogenized(v,h)); }
	std::shared_ptr<Node> CosOperator::Homogenized(VariableGroup const& v, std::shared_ptr<Variable> const& h)    const { return CosOperator::Make(operand_->Homogenized(v,h)); }
	std::shared_ptr<Node> ArcCosOperator::Homogenized(VariableGroup const& v, std::shared_ptr<Variable> const& h) const { return ArcCosOperator::Make(operand_->Homogenized(v,h)); }
	std::shared_ptr<Node> TanOperator::Homogenized(VariableGroup const& v, std::shared_ptr<Variable> const& h)    const { return TanOperator::Make(operand_->Homogenized(v,h)); }
	std::shared_ptr<Node> ArcTanOperator::Homogenized(VariableGroup const& v, std::shared_ptr<Variable> const& h) const { return ArcTanOperator::Make(operand_->Homogenized(v,h)); }
} // re: namespace node
} // re: bertini namespace
