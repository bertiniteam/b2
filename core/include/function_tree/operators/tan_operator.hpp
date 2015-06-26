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
//  Created by Collins, James B. on 6/5/15.
//
//
// tan_operator.hpp:  Declares the class  TanOperator.


#ifndef Function_Tree_Test_tan_operator_hpp
#define Function_Tree_Test_tan_operator_hpp

#include "function_tree/Node.hpp"
#include "function_tree/operators/unary_operator.hpp"
#include "function_tree/operators/mult_operator.hpp"
#include "function_tree/operators/cos_operator.hpp"
#include "function_tree/operators/power_operator.hpp"


namespace bertini {
	
	// Node -> UnaryOperator -> TanOperator
	// Description: This class represents the tangent function.  FreshEval method
	// is defined for tangent and takes the tangent of the child node.
	class TanOperator : public  virtual UnaryOperator
	{
	public:
		
		TanOperator(){}
		
		TanOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		// These do nothing for a constant
		std::string PrintNode() override
		{
			return "-" + child_->PrintNode();
		}
		
		
		virtual void print(std::ostream & target) const override
		{
			target << "tan(";
			child_->print(target);
			target << ")";
		}
        
        
        /**
         Differentiates the tangent function.
         */
        virtual std::shared_ptr<Node> Differentiate() override
        {
            auto ret_mult = std::make_shared<MultOperator>();
            ret_mult->AddChild(child_->Differentiate());
            ret_mult->AddChild(std::make_shared<CosOperator>(child_),false);
            ret_mult->AddChild(std::make_shared<CosOperator>(child_),false);
            return ret_mult;
        }

		
		virtual ~TanOperator() = default;
		
	protected:
		// Specific implementation of FreshEval for negate.
//		dbl FreshEval(dbl) override
//		{
//			return tan(child_->Eval<dbl>());
//		}
//		
//		mpfr FreshEval(mpfr) override
//		{
//			return tan(child_->Eval<mpfr>());
//		}
        
        // Specific implementation of FreshEval for negate.
        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            return tan(child_->Eval<dbl>(diff_variable));
        }
        
        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            return tan(child_->Eval<mpfr>(diff_variable));
        }

	};
	
	
} // re: namespace bertini



namespace {
	// begin the overload of operators
	
	
	inline std::shared_ptr<bertini::Node> tan(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::TanOperator>(N);
	}
	
}

#endif
