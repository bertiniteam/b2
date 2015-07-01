#include "function_tree/operators/intpow_operator.hpp"
#include "function_tree/operators/power_operator.hpp"




namespace bertini {

/**
 Virtual polymorphic method for printing to an arbitrary stream.
 */
 void IntegerPowerOperator::print(std::ostream & target) const
 {
 	target << "(";
 	            child_->print(target);
 	            target << "^" << exponent() << ")";
}


std::shared_ptr<Node> IntegerPowerOperator::Differentiate()
{

	if (exponent_==0.0)
		return std::make_shared<Number>(0.0);
	else if (exponent_==1.0)
		return child_->Differentiate();
	else if (exponent_==2){
		auto M = std::make_shared<MultOperator>(std::make_shared<Number>(2.0), child_);
		M->AddChild(child_->Differentiate());
		return M;
	}
	else{
		auto M = std::make_shared<MultOperator>(std::make_shared<Number>(exponent_), 
		                                        std::make_shared<IntegerPowerOperator>(child_, exponent_-1)
		                                        );
		M->AddChild(child_->Differentiate());
		return M;
	}
}

int IntegerPowerOperator::Degree(std::shared_ptr<Variable> const& v) const
{
	auto base_deg = child_->Degree(v);
	if (base_deg<0)
		return base_deg;
	else
		return exponent_*base_deg;

}


}
