#include "function_tree/operators/trig.hpp"


namespace bertini {
	
	
	std::shared_ptr<Node> SinOperator::Differentiate()
	{
		auto ret_mult = std::make_shared<MultOperator>();
		auto cos_op = std::make_shared<CosOperator>(child_);
		ret_mult->AddChild(cos_op);
		ret_mult->AddChild(child_->Differentiate());
		return ret_mult;
	}
	
	void CosOperator::print(std::ostream & target) const
	{
		target << "cos(";
		child_->print(target);
		target << ")";
	}
	
	std::shared_ptr<Node> CosOperator::Differentiate()
	{
		auto ret_mult = std::make_shared<MultOperator>();
		std::shared_ptr<Node> sin_op = std::make_shared<SinOperator>(child_);
		ret_mult->AddChild(sin_op);
		ret_mult->AddChild(child_->Differentiate());
		return std::make_shared<NegateOperator>(ret_mult);
	}
	
	int CosOperator::Degree(std::shared_ptr<Variable> const& v) const
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
	

	std::shared_ptr<Node> ArcSinOperator::Differentiate()
	{
		return child_->Differentiate()/sqrt(1-pow(child_,2));
	}

	std::shared_ptr<Node> ArcCosOperator::Differentiate()
	{
		return -child_->Differentiate()/sqrt(1-pow(child_,2));
	}


	std::shared_ptr<Node> ArcTanOperator::Differentiate()
	{
		auto ret_mult = std::make_shared<MultOperator>();
		ret_mult->AddChild(child_->Differentiate());
		ret_mult->AddChild(1 + pow(child_,2), false);
		return ret_mult;
	}

	
} // re: bertini namespace
