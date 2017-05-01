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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame
// Jeb Collins, West Texas A&M


#include "function_tree/operators/trig.hpp"


namespace bertini {
namespace node{	

	int TrigOperator::Degree(std::shared_ptr<Variable> const& v) const
	{
		if (child_->Degree(v)==0)
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
		child_->print(target);
		target << ")";
	}



	std::shared_ptr<Node> SinOperator::Differentiate() const
	{
		return cos(child_) * child_->Differentiate();
	}
	
	std::shared_ptr<Node> SinOperator::Derivative(std::shared_ptr<Variable> const& v) const
	{
		return cos(child_) * child_->Derivative(v);
	}

	// Specific implementation of FreshEval for negate.
	dbl SinOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return sin(child_->Eval<dbl>(diff_variable));
	}
	
	mpfr SinOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return sin(child_->Eval<mpfr>(diff_variable));
	}
	
	
	void SinOperator::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<mpfr>(evaluation_value, diff_variable);
		evaluation_value = sin(evaluation_value);
	}
	
	void SinOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = sin(evaluation_value);
	}

	void ArcSinOperator::print(std::ostream & target) const
	{
		target << "asin(";
		child_->print(target);
		target << ")";
	}


	std::shared_ptr<Node> ArcSinOperator::Differentiate() const
	{
		return child_->Differentiate()/sqrt(1-pow(child_,2));
	}

	std::shared_ptr<Node> ArcSinOperator::Derivative(std::shared_ptr<Variable> const& v) const
	{
		return child_->Derivative(v)/sqrt(1-pow(child_,2));
	}


	// Specific implementation of FreshEval for negate.
	dbl ArcSinOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return asin(child_->Eval<dbl>(diff_variable));
	}
	
	mpfr ArcSinOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return asin(child_->Eval<mpfr>(diff_variable));
	}

	
	void ArcSinOperator::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<mpfr>(evaluation_value, diff_variable);
		evaluation_value = asin(evaluation_value);
	}
	
	void ArcSinOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = asin(evaluation_value);
	}

	void CosOperator::print(std::ostream & target) const
	{
		target << "cos(";
		child_->print(target);
		target << ")";
	}
	
	std::shared_ptr<Node> CosOperator::Differentiate() const
	{
		return -sin(child_) * child_->Differentiate();
	}
	
	std::shared_ptr<Node> CosOperator::Derivative(std::shared_ptr<Variable> const& v) const
	{
		return -sin(child_) * child_->Derivative(v);
	}

	dbl CosOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return cos(child_->Eval<dbl>(diff_variable));
	}
	
	mpfr CosOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return cos(child_->Eval<mpfr>(diff_variable));
	}
	
	
	void CosOperator::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<mpfr>(evaluation_value, diff_variable);
		evaluation_value = cos(evaluation_value);
	}
	
	void CosOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = cos(evaluation_value);
	}




	void ArcCosOperator::print(std::ostream & target) const
	{
		target << "acos(";
		child_->print(target);
		target << ")";
	}




	std::shared_ptr<Node> ArcCosOperator::Differentiate() const
	{
		return -child_->Differentiate()/sqrt(1-pow(child_,2));
	}

	std::shared_ptr<Node> ArcCosOperator::Derivative(std::shared_ptr<Variable> const& v) const
	{
		return -child_->Derivative(v)/sqrt(1-pow(child_,2));
	}

	// Specific implementation of FreshEval for negate.
	dbl ArcCosOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return acos(child_->Eval<dbl>(diff_variable));
	}
	
	mpfr ArcCosOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return acos(child_->Eval<mpfr>(diff_variable));
	}
	
	
	void ArcCosOperator::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<mpfr>(evaluation_value, diff_variable);
		evaluation_value = acos(evaluation_value);
	}
	
	void ArcCosOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = acos(evaluation_value);
	}



	void TanOperator::print(std::ostream & target) const
	{
		target << "tan(";
		child_->print(target);
		target << ")";
	}


	std::shared_ptr<Node> TanOperator::Differentiate() const
	{
		return child_->Differentiate() /  pow(cos(child_),2);
	}

	std::shared_ptr<Node> TanOperator::Derivative(std::shared_ptr<Variable> const& v) const
	{
		return child_->Derivative(v) /  pow(cos(child_),2);
	}

	dbl TanOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return tan(child_->Eval<dbl>(diff_variable));
	}
	
	mpfr TanOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return tan(child_->Eval<mpfr>(diff_variable));
	}
	
	
	void TanOperator::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<mpfr>(evaluation_value, diff_variable);
		evaluation_value = tan(evaluation_value);
	}
	
	void TanOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = tan(evaluation_value);
	}




	void ArcTanOperator::print(std::ostream & target) const
	{
		target << "atan(";
		child_->print(target);
		target << ")";
	}
		


	std::shared_ptr<Node> ArcTanOperator::Differentiate() const
	{
		return child_->Differentiate() / (1 + pow(child_,2));
	}

	std::shared_ptr<Node> ArcTanOperator::Derivative(std::shared_ptr<Variable> const& v) const
	{
		return child_->Derivative(v) / (1 + pow(child_,2));
	}

	dbl ArcTanOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
	{
		return atan(child_->Eval<dbl>(diff_variable));
	}
	
	mpfr ArcTanOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
	{
		return atan(child_->Eval<mpfr>(diff_variable));
	}
	
	
	void ArcTanOperator::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<mpfr>(evaluation_value, diff_variable);
		evaluation_value = atan(evaluation_value);
	}
	
	void ArcTanOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
	{
		child_->EvalInPlace<dbl>(evaluation_value, diff_variable);
		evaluation_value = atan(evaluation_value);
	}
} // re: namespace node	
} // re: bertini namespace
