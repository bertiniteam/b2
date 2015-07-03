#include "function_tree/operators/arithmetic.hpp"


namespace bertini{
	
	using ::pow;
	
	///////////////////////
	//
	//  Sum Operator definitions
	//
	//////////////////////
	
	void SumOperator::print(std::ostream & target) const
	{
		target << "(";
		for (auto iter = children_.begin(); iter!= children_.end(); iter++) {
			if (iter==children_.begin()) {
				// on the first iteration, no need to put a + if a +
				if ( !(*(children_sign_.begin()+(iter-children_.begin()))) )
					target << "-";
			}
			else
			{
				if ( !(*(children_sign_.begin()+(iter-children_.begin()))) )
					target << "-";
				else
					target << "+";
			}
			(*iter)->print(target);
			
		}
		target << ")";
	}
	
	std::shared_ptr<Node> SumOperator::Differentiate()
	{
		unsigned int counter = 0;
		std::shared_ptr<SumOperator> ret_sum = std::make_shared<SumOperator>();
		for (int ii = 0; ii < children_.size(); ++ii)
		{
			auto converted = std::dynamic_pointer_cast<Number>(children_[ii]);
			if (converted)
				continue;
			
			auto temp_node = children_[ii]->Differentiate();
			converted = std::dynamic_pointer_cast<Number>(temp_node);
			if (converted)
				if (converted->Eval<dbl>()==dbl(0.0))
					continue;
			
			ret_sum->AddChild(temp_node,children_sign_[ii]);
			counter++;
		}
		
		if (counter>0)
			return ret_sum;
		else
			return std::make_shared<Number>(0.0);
	}
	
	
	
	
	int SumOperator::Degree(std::shared_ptr<Variable> const& v) const
	{
		int deg = 0;
		
		for (auto iter: children_)
		{
			auto curr_deg = iter->Degree(v);
			if (curr_deg<0)
				return curr_deg;
			
			deg = std::max(deg, curr_deg);
		}
		return deg;
	}
	
	int SumOperator::Degree(VariableGroup const& vars) const 
	{
		auto deg = 0;
		for (auto iter = children_.begin(); iter!=children_.end(); iter++)
		{
			auto term_degree = (*iter)->Degree(vars);
			if (term_degree<0)
				return term_degree;

			deg = std::max(deg, term_degree);
			
		}
		
		return deg;
	}


	std::vector<int> SumOperator::MultiDegree(VariableGroup const& vars) const
	{
		std::vector<int> deg(vars.size(),0);
		for (auto iter : children_)
		{
			auto term_deg = iter->MultiDegree(vars);

			for (auto iter = term_deg.begin(); iter!= term_deg.end(); ++iter)
			{
				*(deg.begin()+(iter-term_deg.begin())) = std::max(*(deg.begin()+(iter-term_deg.begin())), *iter);
			}
		}
		return deg;
	}

	void SumOperator::Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar)
	{
		
		
		
		// first homogenize each summand.
		for (auto iter: children_)
		{
			iter->Homogenize(vars, homvar);
		}
		
		// then, homogenize this sum.
		
		// compute the highest degree among all summands.
		int maxdegree = 0;
		std::vector<int> term_degrees;
		// first homogenize each summand.
		for (auto iter: children_)
		{
			auto local_degree = iter->Degree(vars);
			if (local_degree<0)
				throw std::runtime_error("asking for homogenization on non-polynomial node");
			// TODO: this throw would leave the tree in a partially homogenized state.  this is scary.
			
			term_degrees.push_back(local_degree);
			maxdegree = std::max(maxdegree, local_degree);
		}
		
		for (auto iter = children_.begin(); iter!=children_.end(); iter++)
		{
			auto degree_deficiency = maxdegree - *(term_degrees.begin() + (iter-children_.begin()));
			if ( degree_deficiency > 0)
			{

				// hold the child temporarily.
				if (degree_deficiency==1)
				{
					std::shared_ptr<Node> M = std::make_shared<MultOperator>(homvar,std::dynamic_pointer_cast<Node>(*iter));
					swap(*iter,M);
				}
				else{
					std::shared_ptr<Node> P = std::make_shared<IntegerPowerOperator>(std::dynamic_pointer_cast<Node>(homvar),degree_deficiency);
					std::shared_ptr<Node> M = std::make_shared<MultOperator>(P,std::dynamic_pointer_cast<Node>(*iter));
					swap(*iter,M);
				}
				

			}
		}
		
	}


	bool SumOperator::IsHomogeneous(std::shared_ptr<Variable> const& v) const
	{
		
		for (auto iter : children_)
		{
			if (!iter->IsHomogeneous(v))
				return false;
		}

		// the only hope this has of being homogeneous, is that each factor is homogeneous
		int deg;

		deg = (*(children_.begin()))->Degree(v) ;
	
		if (deg < 0) 
			return false;

		for (auto iter = children_.begin()+1; iter!= children_.end(); iter++)
		{
			auto local_degree = (*iter)->Degree(v);
			if (local_degree!=deg)
				return false;
		}

		return true;
	}

	bool SumOperator::IsHomogeneous(VariableGroup const& v) const
	{
		
		for (auto iter : children_)
		{
			if (!iter->IsHomogeneous(v))
				return false;
		}

		// the only hope this has of being homogeneous, is that each factor is homogeneous
		int deg;

		deg = (*(children_.begin()))->Degree(v) ;
	
		if (deg < 0) 
			return false;

		for (auto iter = children_.begin()+1; iter!= children_.end(); iter++)
		{
			auto local_degree = (*iter)->Degree(v);
			if (local_degree!=deg)
				return false;
		}

		return true;
	}

	

	
	dbl SumOperator::FreshEval(dbl, std::shared_ptr<Variable> diff_variable)
	{
		dbl retval{0};
		for(int ii = 0; ii < children_.size(); ++ii)
		{
			if(children_sign_[ii])
			{
				retval += children_[ii]->Eval<dbl>(diff_variable);
			}
			else
			{
				retval -= children_[ii]->Eval<dbl>(diff_variable);
			}
		}
		
		return retval;
	}
	
	mpfr SumOperator::FreshEval(mpfr, std::shared_ptr<Variable> diff_variable)
	{
		mpfr retval{0};
		for(int ii = 0; ii < children_.size(); ++ii)
		{
			if(children_sign_[ii])
			{
				retval += children_[ii]->Eval<mpfr>(diff_variable);
			}
			else
			{
				retval -= children_[ii]->Eval<mpfr>(diff_variable);
			}
		}
		
		return retval;
	}
	
	
	















	//////////////////////
	//
	//  Negate operator definitions
	//
	////////////////////////
	
	
	void NegateOperator::print(std::ostream & target) const
	{
		target << "-(";
		child_->print(target);
		target << ")";
	}
	
	std::shared_ptr<Node> NegateOperator::Differentiate()
	{
		return std::make_shared<NegateOperator>(child_->Differentiate());
	}
	
	dbl NegateOperator::FreshEval(dbl, std::shared_ptr<Variable> diff_variable)
	{
		return -(child_->Eval<dbl>(diff_variable));
	}
	
	mpfr NegateOperator::FreshEval(mpfr, std::shared_ptr<Variable> diff_variable)
	{
		return -child_->Eval<mpfr>(diff_variable);
	}
	
	
	
	
















	
	
	
	///////////////////////
	//
	//  Mult Operator definitions
	//
	//////////////////////
	
	void MultOperator::print(std::ostream & target) const
	{
		target << "(";
		for (auto iter = children_.begin(); iter!= children_.end(); iter++) {
			if (iter==children_.begin())
				if (! *children_mult_or_div_.begin()  )
					target << "1/";
			(*iter)->print(target);
			if (iter!=(children_.end()-1)){
				if (*(children_mult_or_div_.begin() + (iter-children_.begin())+1)) { // TODO i think this +1 is wrong... dab
					target << "*";
				}
				else{
					target << "/";
				}
				
			}
		}
		target << ")";
	}
	
	
	std::shared_ptr<Node> MultOperator::Differentiate()
	{
		std::shared_ptr<SumOperator> ret_sum = std::make_shared<SumOperator>();
		
		unsigned term_counter {0};
		
		for (int ii = 0; ii < children_.size(); ++ii)
		{
			
			auto local_derivative = children_[ii]->Differentiate();
			
			auto is_it_a_number = std::dynamic_pointer_cast<Number>(local_derivative);
			if (is_it_a_number)
				if (is_it_a_number->Eval<dbl>()==dbl(0.0))
					continue;
			
			
			
			auto term_ii = std::make_shared<MultOperator>();
			for (int jj = 0; jj < children_.size(); ++jj)
			{
				if(jj != ii)
					term_ii->AddChild(children_[jj],children_mult_or_div_[jj]);
			}
			
			
			// if the derivative of the term under consideration is equal to 1, no need to go on.  already have the product we need.
			if (is_it_a_number)
				if (is_it_a_number->Eval<dbl>()==dbl(1.0))
					continue;
			
			
			// the first branch is for division, and the quotient rule.
			if ( !(children_mult_or_div_[ii]) )
			{
				// divide by the
				term_ii->AddChild(children_[ii],false);
				term_ii->AddChild(children_[ii],false);
				term_ii->AddChild(local_derivative);
				ret_sum->AddChild(term_ii,false);  // false for subtract here.
			}
			// the second branch is for multiplication, and the product rule
			else
			{
				term_ii->AddChild(local_derivative,true);// true for multiply, not divide
				ret_sum->AddChild(term_ii); //  no truth second argument, because is add.
			}
			
			term_counter++;
			
		} // re: for ii
		
		return ret_sum;
	}
	
	
	
	
	int MultOperator::Degree(std::shared_ptr<Variable> const& v) const
	{
		int deg = 0;
		for (auto iter = children_.begin(); iter!= children_.end(); iter++)
		{
			// if the child node is a differential coming from another variable, then this degree should be 0, end of story.
			
			auto factor_deg = (*iter)->Degree(v);
			
			auto is_it_a_differential = std::dynamic_pointer_cast<Differential>(*iter);
			if (is_it_a_differential)
				if (is_it_a_differential->GetVariable()!=v)
				{
					// std::cout << *this << " deg " << 0 << "\n";
					return 0;
				}
			
			
			
			
			
			if (factor_deg<0)
				return factor_deg;
			else if (factor_deg!=0 && !*(children_mult_or_div_.begin() + (iter-children_.begin()) ) )
				return -1;
			else
				deg+=factor_deg;
		}
		// std::cout << *this << " deg " << deg << "\n";
		return deg;
	}
	

	int MultOperator::Degree(VariableGroup const& vars) const 
	{
		
		auto deg = 0;

		for (auto iter = children_.begin(); iter!=children_.end(); iter++)
		{
			auto factor_deg = (*iter)->Degree(vars);  

			if (factor_deg<0)
				return factor_deg;
			else if (factor_deg!=0 && !*(children_mult_or_div_.begin() + (iter-children_.begin()) ) )
				return -1;
			else
				deg+=factor_deg;
		}

		
		return deg;
	}

	std::vector<int> MultOperator::MultiDegree(VariableGroup const& vars) const
	{
		std::vector<int> deg(vars.size(),0);
		for (auto iter : children_)
		{
			auto term_deg = iter->MultiDegree(vars);

			for (auto iter = term_deg.begin(); iter!= term_deg.end(); ++iter)
			{
				*(deg.begin()+(iter-term_deg.begin())) += *iter;
			}
		}
		return deg;
	}



	void MultOperator::Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar)
	{
		for (auto iter: children_)
		{
			iter->Homogenize(vars, homvar);
		}
	}


	bool MultOperator::IsHomogeneous(std::shared_ptr<Variable> const& v) const
	{
		// the only hope this has of being homogeneous, is that each factor is homogeneous
		for (auto iter : children_)
		{
			if (! iter->IsHomogeneous(v))
			{
				return false;
			}
		}
		return true;
	}
	
	bool MultOperator::IsHomogeneous(VariableGroup const& v) const
	{
		// the only hope this has of being homogeneous, is that each factor is homogeneous
		for (auto iter : children_)
		{
			if (! iter->IsHomogeneous(v))
			{
				return false;
			}
		}
		return true;
	}

	dbl MultOperator::FreshEval(dbl, std::shared_ptr<Variable> diff_variable)
	{
		dbl retval{1};
		for(int ii = 0; ii < children_.size(); ++ii)
		{
			if(children_mult_or_div_[ii])
			{
				retval *= children_[ii]->Eval<dbl>(diff_variable);
			}
			else
			{
				retval /= children_[ii]->Eval<dbl>(diff_variable);
			}
		}
		
		return retval;
	}
	
	mpfr MultOperator::FreshEval(mpfr, std::shared_ptr<Variable> diff_variable)
	{
		mpfr retval{1};
		for(int ii = 0; ii < children_.size(); ++ii)
		{
			if(children_mult_or_div_[ii])
			{
				retval *= children_[ii]->Eval<mpfr>(diff_variable);
			}
			else
			{
				retval /= children_[ii]->Eval<mpfr>(diff_variable);
			}
		}
		
		return retval;
	}
	
	
	////////////
	//
	//  Power Operator definitions
	//
	/////////////////
	
	void PowerOperator::Reset()
	{
		Node::ResetStoredValues();
		base_->Reset();
		exponent_->Reset();
	}
	
	void PowerOperator::print(std::ostream & target) const
	{
		target << "(" << *base_ << ")^(" << *exponent_ << ")";
	}
	
	
	std::shared_ptr<Node> PowerOperator::Differentiate()
	{
		auto ret_mult = std::make_shared<MultOperator>();
		auto exp_minus_one = std::make_shared<SumOperator>(exponent_, true, std::make_shared<Number>("1.0"),false);
		ret_mult->AddChild(base_->Differentiate());
		ret_mult->AddChild(exponent_);
		ret_mult->AddChild(std::make_shared<PowerOperator>(base_, exp_minus_one));
		return ret_mult;
	}
	
	
	int PowerOperator::Degree(std::shared_ptr<Variable> const& v) const
	{
		
		auto base_deg = base_->Degree(v);
		auto exp_deg = exponent_->Degree(v);
		
		if (exp_deg==0)
		{
			auto exp_val = exponent_->Eval<dbl>();
			bool exp_is_int = false;
			
			if (fabs(imag(exp_val))< 10*std::numeric_limits<double>::epsilon()) // so a real thresholding step
				if (fabs(real(exp_val) - std::round(real(exp_val))) < 10*std::numeric_limits<double>::epsilon()) // then check that the real part is close to an integer
					exp_is_int = true;
			
			if (exp_is_int)
			{
				if (abs(exp_val-dbl(0.0))< 10*std::numeric_limits<double>::epsilon())
					return 0;
				else if (real(exp_val)<0)
					return -1;
				else  // positive integer.
				{
					if (base_deg<0)
						return -1;
					else
						return base_deg*std::round(real(exp_val));
				}
				
			}
			else
			{
				if (base_deg==0)
					return 0;
				else
					return -1;
			}
			
		}
		else
		{
			// there may be an edge case here where the base is the numbers 0 or 1.  but that would be stupid, wouldn't it.
			return -1;
		}
	}
	
	int PowerOperator::Degree(VariableGroup const& vars) const
	{
		auto multideg = MultiDegree(vars);
		auto deg = 0;
		std::for_each(multideg.begin(),multideg.end(),[&](int n){
						if (n < 0)
							deg = -1;
						else
                        	deg += n;
 						});
		return deg;
	}


	std::vector<int> PowerOperator::MultiDegree(VariableGroup const& vars) const
	{
		std::vector<int> deg(vars.size(),0);
		for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
		{
			*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
		}
		return deg;
	}



	void PowerOperator::Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar)
	{
		if (exponent_->Degree(vars)==0)
		{
			base_->Homogenize(vars, homvar);
		}
		else{
			throw std::runtime_error("asking for homogenization on non-polynomial node");
			//TODO: this will leave the system in a broken state, partially homogenized...
		}
	}

	bool PowerOperator::IsHomogeneous(std::shared_ptr<Variable> const& v) const
	{
		// the only hope this has of being homogeneous, is that the degree of the exponent is 0 (it's constant), and that it's an integer
		if (exponent_->Degree(v)==0)
		{
			auto exp_val = exponent_->Eval<dbl>();
			if (fabs(imag(exp_val)) < 10*std::numeric_limits<double>::epsilon())
				if (fabs(std::round(real(exp_val)) - real(exp_val)) < 10*std::numeric_limits<double>::epsilon())
					if (real(exp_val) >=0 )
						return base_->IsHomogeneous(v);
		}
		return false;
	}


	bool PowerOperator::IsHomogeneous(VariableGroup const& v) const
	{
		// the only hope this has of being homogeneous, is that the degree of the exponent is 0 (it's constant), and that it's an integer
		if (exponent_->Degree(v)==0)
		{
			auto exp_val = exponent_->Eval<dbl>();
			if (fabs(imag(exp_val)) < 10*std::numeric_limits<double>::epsilon())
				if (fabs(std::round(real(exp_val)) - real(exp_val)) < 10*std::numeric_limits<double>::epsilon())
					if (real(exp_val) >=0 )
						return base_->IsHomogeneous(v);
		}
		return false;
	}
	
	dbl PowerOperator::FreshEval(dbl, std::shared_ptr<Variable> diff_variable)
	{
		return std::pow( base_->Eval<dbl>(diff_variable), exponent_->Eval<dbl>());
	}
	
	mpfr PowerOperator::FreshEval(mpfr, std::shared_ptr<Variable> diff_variable)
	{
		return pow( base_->Eval<mpfr>(diff_variable), exponent_->Eval<mpfr>());
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/////////////////
	//
	//  IntegerPowerOperator definitions
	//
	////////////////////
	
	void IntegerPowerOperator::print(std::ostream & target) const
	{
		target << "(";
		child_->print(target);
		target << "^" << exponent() << ")";
	}
	
	
	std::shared_ptr<Node> IntegerPowerOperator::Differentiate()
	{
		
		if (exponent_==0)
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
													std::make_shared<IntegerPowerOperator>(child_, exponent_-1) );
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
	

	
	
	
	
	
	
	
	
	
	
	//////////////
	//
	//  Square Root Operator definitions
	//
	/////////////////
	
	void SqrtOperator::print(std::ostream & target) const
	{
		target << "sqrt(";
		child_->print(target);
		target << ")";
	}
	
	
	
	std::shared_ptr<Node> SqrtOperator::Differentiate()
	{
		auto ret_mult = std::make_shared<MultOperator>();
		ret_mult->AddChild(std::make_shared<PowerOperator>(child_, std::make_shared<Number>(-0.5)));
		ret_mult->AddChild(child_->Differentiate());
		ret_mult->AddChild(std::make_shared<Number>(0.5));
		return ret_mult;
	}
	
	int SqrtOperator::Degree(std::shared_ptr<Variable> const& v) const
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
	

	dbl SqrtOperator::FreshEval(dbl, std::shared_ptr<Variable> diff_variable)
	{
		return sqrt(child_->Eval<dbl>(diff_variable));
	}
	
	mpfr SqrtOperator::FreshEval(mpfr, std::shared_ptr<Variable> diff_variable)
	{
		return sqrt(child_->Eval<mpfr>(diff_variable));
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	///////////////
	//
	//  ExpOperator definitions
	//
	//////////////
	
	void ExpOperator::print(std::ostream & target) const
	{
		target << "exp(";
		child_->print(target);
		target << ")";
	}
	
	std::shared_ptr<Node> ExpOperator::Differentiate()
	{
		auto ret_mult = std::make_shared<MultOperator>();
		ret_mult->AddChild(std::make_shared<ExpOperator>(child_));
		ret_mult->AddChild(child_->Differentiate());
		return ret_mult;
	}
	
	int ExpOperator::Degree(std::shared_ptr<Variable> const& v) const
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
	
	dbl ExpOperator::FreshEval(dbl, std::shared_ptr<Variable> diff_variable)
	{
		return exp(child_->Eval<dbl>(diff_variable));
	}
	
	mpfr ExpOperator::FreshEval(mpfr, std::shared_ptr<Variable> diff_variable)
	{
		return exp(child_->Eval<mpfr>(diff_variable));
	}
	
}
