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
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
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


	unsigned SinOperator::EliminateZeros()
	{
		return 0;
	}
	unsigned SinOperator::EliminateOnes()
	{
		return 0;
	}

	void SinOperator::print(std::ostream & target) const
	{
		target << "sin(";
		operand_->print(target);
		target << ")";
	}

	
	std::shared_ptr<Node> SinOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return cos(operand_) * operand_->Differentiate(v);
	}

	// Specific implementation of FreshEval for negate.
	dbl SinOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return sin(operand_->Eval<dbl>(diff_variable));
	}
	
	mpfr_complex SinOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return sin(operand_->Eval<mpfr_complex>(diff_variable));
	}
	
	
	void SinOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
		evaluation_value = sin(evaluation_value);
	}
	
	void SinOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = sin(evaluation_value);
	}

	unsigned ArcSinOperator::EliminateZeros()
	{
		return 0;
	}
	unsigned ArcSinOperator::EliminateOnes()
	{
		return 0;
	}

	void ArcSinOperator::print(std::ostream & target) const
	{
		target << "asin(";
		operand_->print(target);
		target << ")";
	}

	std::shared_ptr<Node> ArcSinOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return operand_->Differentiate(v)/sqrt(1-pow(operand_,2));
	}


	// Specific implementation of FreshEval for negate.
	dbl ArcSinOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return asin(operand_->Eval<dbl>(diff_variable));
	}
	
	mpfr_complex ArcSinOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return asin(operand_->Eval<mpfr_complex>(diff_variable));
	}

	
	void ArcSinOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
		evaluation_value = asin(evaluation_value);
	}
	
	void ArcSinOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = asin(evaluation_value);
	}

	unsigned CosOperator::EliminateZeros()
	{
		return 0;
	}
	unsigned CosOperator::EliminateOnes()
	{
		return 0;
	}

	void CosOperator::print(std::ostream & target) const
	{
		target << "cos(";
		operand_->print(target);
		target << ")";
	}
	
	std::shared_ptr<Node> CosOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return -sin(operand_) * operand_->Differentiate(v);
	}

	dbl CosOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return cos(operand_->Eval<dbl>(diff_variable));
	}
	
	mpfr_complex CosOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return cos(operand_->Eval<mpfr_complex>(diff_variable));
	}
	
	
	void CosOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
		evaluation_value = cos(evaluation_value);
	}
	
	void CosOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = cos(evaluation_value);
	}


	unsigned ArcCosOperator::EliminateZeros()
	{
		return 0;
	}
	unsigned ArcCosOperator::EliminateOnes()
	{
		return 0;
	}

	void ArcCosOperator::print(std::ostream & target) const
	{
		target << "acos(";
		operand_->print(target);
		target << ")";
	}


	std::shared_ptr<Node> ArcCosOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return -operand_->Differentiate(v)/sqrt(1-pow(operand_,2));
	}

	// Specific implementation of FreshEval for negate.
	dbl ArcCosOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return acos(operand_->Eval<dbl>(diff_variable));
	}
	
	mpfr_complex ArcCosOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return acos(operand_->Eval<mpfr_complex>(diff_variable));
	}
	
	
	void ArcCosOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
		evaluation_value = acos(evaluation_value);
	}
	
	void ArcCosOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = acos(evaluation_value);
	}


	unsigned TanOperator::EliminateZeros()
	{
		return 0;
	}
	unsigned TanOperator::EliminateOnes()
	{
		return 0;
	}

	void TanOperator::print(std::ostream & target) const
	{
		target << "tan(";
		operand_->print(target);
		target << ")";
	}

	std::shared_ptr<Node> TanOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return operand_->Differentiate(v) /  pow(cos(operand_),2);
	}

	dbl TanOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return tan(operand_->Eval<dbl>(diff_variable));
	}
	
	mpfr_complex TanOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return tan(operand_->Eval<mpfr_complex>(diff_variable));
	}
	
	
	void TanOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
		evaluation_value = tan(evaluation_value);
	}
	
	void TanOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = tan(evaluation_value);
	}


	unsigned ArcTanOperator::EliminateZeros()
	{
		return 0;
	}
	unsigned ArcTanOperator::EliminateOnes()
	{
		return 0;
	}

	void ArcTanOperator::print(std::ostream & target) const
	{
		target << "atan(";
		operand_->print(target);
		target << ")";
	}
		

	std::shared_ptr<Node> ArcTanOperator::Differentiate(std::shared_ptr<Variable> const& v) const
	{
		return operand_->Differentiate(v) / (1 + pow(operand_,2));
	}

	dbl ArcTanOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return atan(operand_->Eval<dbl>(diff_variable));
	}
	
	mpfr_complex ArcTanOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return atan(operand_->Eval<mpfr_complex>(diff_variable));
	}
	
	
	void ArcTanOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
		evaluation_value = atan(evaluation_value);
	}
	
	void ArcTanOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = atan(evaluation_value);
	}
} // re: namespace node	
} // re: bertini namespace
