//This file is part of Bertini 2.0.
//
//exp_operator.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//exp_operator.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with exp_operator.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// exp_operator.hpp:  Declares the class  ExpOperator.



#ifndef Function_Tree_Test_exp_operator_hpp
#define Function_Tree_Test_exp_operator_hpp

#include "function_tree/Node.hpp"
#include "function_tree/operators/unary_operator.hpp"
#include "function_tree/operators/mult_operator.hpp"

namespace bertini {
	
	
	
	// Node -> UnaryOperator -> ExpOperator
	// Description: This class represents the exponential function.  FreshEval method
	// is defined for exponential and takes the exponential of the child node.
	class ExpOperator : public  virtual UnaryOperator
	{
	public:
		
		ExpOperator(){}
		
		ExpOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
	 
		// These do nothing for a constant
		std::string PrintNode() override
		{
			return "-" + child_->PrintNode();
		}
		
		
		virtual void print(std::ostream & target) const override
		{
			target << "exp(";
			child_->print(target);
			target << ")";
		}
        
        
        /**
         Differentiates the exponential function.
         */
        virtual std::shared_ptr<Node> Differentiate() override
        {
            auto ret_mult = std::make_shared<MultOperator>();
            ret_mult->AddChild(std::make_shared<ExpOperator>(child_));
            ret_mult->AddChild(child_->Differentiate());
            return ret_mult;
        }

		
		virtual ~ExpOperator() = default;
		
	protected:
		// Specific implementation of FreshEval for exponentiate.
		dbl FreshEval(dbl) override
		{
			return exp(child_->Eval<dbl>());
		}
		
		mpfr FreshEval(mpfr) override
		{
			return exp(child_->Eval<mpfr>());
		}
        
        // Specific implementation of FreshEval for exponentiate.
        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            return exp(child_->Eval<dbl>(diff_variable));
        }
        
        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            return exp(child_->Eval<mpfr>(diff_variable));
        }

	};

	
} // re: namespace bertini


namespace {
	// begin the overload of operators
	
	
	inline std::shared_ptr<bertini::Node> exp(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::ExpOperator>(N);
	}
}

#endif
