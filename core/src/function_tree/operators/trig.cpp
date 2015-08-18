//This file is part of Bertini 2.0.
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

//  trig.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015


#include "function_tree/operators/trig.hpp"


namespace bertini {
	

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



	std::shared_ptr<Node> SinOperator::Differentiate()
	{
		return cos(child_) * child_->Differentiate();
	}
	



	void ArcSinOperator::print(std::ostream & target) const
	{
		target << "asin(";
		child_->print(target);
		target << ")";
	}


	std::shared_ptr<Node> ArcSinOperator::Differentiate()
	{
		return child_->Differentiate()/sqrt(1-pow(child_,2));
	}




	void CosOperator::print(std::ostream & target) const
	{
		target << "cos(";
		child_->print(target);
		target << ")";
	}
	
	std::shared_ptr<Node> CosOperator::Differentiate()
	{
		return -sin(child_) * child_->Differentiate();
	}
	
	






	void ArcCosOperator::print(std::ostream & target) const
	{
		target << "acos(";
		child_->print(target);
		target << ")";
	}




	std::shared_ptr<Node> ArcCosOperator::Differentiate()
	{
		return -child_->Differentiate()/sqrt(1-pow(child_,2));
	}







	void TanOperator::print(std::ostream & target) const
	{
		target << "tan(";
		child_->print(target);
		target << ")";
	}


	std::shared_ptr<Node> TanOperator::Differentiate()
	{
		return child_->Differentiate() /  pow(cos(child_),2);
	}






	void ArcTanOperator::print(std::ostream & target) const
	{
		target << "atan(";
		child_->print(target);
		target << ")";
	}
		


	std::shared_ptr<Node> ArcTanOperator::Differentiate()
	{
		return child_->Differentiate() / (1 + pow(child_,2));
	}

	
} // re: bertini namespace
