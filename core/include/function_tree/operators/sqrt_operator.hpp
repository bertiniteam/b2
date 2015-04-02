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
// sqrt_operator.hpp:  Declares the class  SqrtOperator.


#ifndef Function_Tree_Test_sqrt_operator_hpp
#define Function_Tree_Test_sqrt_operator_hpp

#include "function_tree/Node.hpp"
#include "function_tree/operators/unary_operator.hpp"


namespace bertini {
    
    // Node -> UnaryOperator -> SqrtOperator
    // Description: This class represents the square root function.  FreshEval method
    // is defined for square root and takes the square root of the child node.
    class SqrtOperator : public  virtual UnaryOperator
    {
    public:
        
        SqrtOperator(){}
        
        SqrtOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
        {};
        
        // These do nothing for a constant
        std::string PrintNode() override
        {
            return "-" + child_->PrintNode();
        }
        
        
        virtual void print(std::ostream & target) const override
        {
        target << "sqrt(";
        child_->print(target);
        target << ")";
        }
        
        virtual ~SqrtOperator() = default;
        
    protected:
        // Specific implementation of FreshEval for negate.
        dbl FreshEval(dbl) override
        {
            return sqrt(child_->Eval<dbl>());
        }
        
        mpfr FreshEval(mpfr) override
        {
            return sqrt(child_->Eval<mpfr>());
        }
        };
        
        
        // begin the overload of operators
        
        
        inline std::shared_ptr<bertini::Node> sqrt(const std::shared_ptr<bertini::Node> & N)
        {
            return std::make_shared<bertini::SqrtOperator>(N);
        }
        
        
        
        
        } // re: namespace bertini
        
        
        
#endif
