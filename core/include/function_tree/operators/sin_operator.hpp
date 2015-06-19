//This file is part of Bertini 2.0.
//
//negate_operator.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//negate_operator.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with negate_operator.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 6/4/15.
//
//
// sin_operator.hpp:  Declares the class  SinOperator.


#ifndef Function_Tree_Test_sin_operator_hpp
#define Function_Tree_Test_sin_operator_hpp

#include "function_tree/node.hpp"
#include "function_tree/operators/unary_operator.hpp"
#include "function_tree/operators/cos_operator.hpp"




namespace bertini {
    
    
    
    // Node -> UnaryOperator -> SinOperator
    // Description: This class represents the sine function.  FreshEval method
    // is defined for sine and takes the sine of the child node.
    class SinOperator : public  virtual UnaryOperator
    {
    public:
        
        SinOperator(){}
        
        SinOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
        {};
        
        // These do nothing for a constant
        std::string PrintNode() override
        {
            return "-" + child_->PrintNode();
        }
        
        
        virtual void print(std::ostream & target) const override
        {
            target << "sin(";
            child_->print(target);
            target << ")";
        }
        
        
        
        /**
         Differentiates the sine function.
         */
        virtual std::shared_ptr<Node> Differentiate() override
        {
            auto ret_mult = std::make_shared<MultOperator>();
            auto cos_op = std::make_shared<CosOperator>(child_);
            ret_mult->AddChild(cos_op);
            ret_mult->AddChild(child_->Differentiate());
            return ret_mult;
        }

        
        virtual ~SinOperator() = default;
        
    protected:
        // Specific implementation of FreshEval for negate.
//        dbl FreshEval(dbl) override
//        {
//            return sin(child_->Eval<dbl>());
//        }
//        
//        mpfr FreshEval(mpfr) override
//        {
//            return sin(child_->Eval<mpfr>());
//        }
        
        // Specific implementation of FreshEval for negate.
        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            return sin(child_->Eval<dbl>(diff_variable));
        }
        
        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            return sin(child_->Eval<mpfr>(diff_variable));
        }

    };
    

    
} // re: namespace bertini




namespace {
	// begin the overload of operators
	
	
	inline std::shared_ptr<bertini::Node> sin(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::SinOperator>(N);
	}
}



#endif
