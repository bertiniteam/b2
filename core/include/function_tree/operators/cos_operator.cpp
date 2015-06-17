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
// cos_operator.cpp:  Defines differentiation for CosOperator.


#include "function_tree/operators/cos_operator.hpp"
#include "function_tree/operators/sin_operator.hpp"



namespace bertini {

std::shared_ptr<Node> CosOperator::Differentiate()
{
    auto ret_mult = std::make_shared<MultOperator>();
    std::shared_ptr<Node> sin_op = std::make_shared<SinOperator>(child_);
    ret_mult->AddChild(sin_op);
    ret_mult->AddChild(child_->Differentiate());
    return std::make_shared<NegateOperator>(ret_mult);
}



} // re: bertini namespace