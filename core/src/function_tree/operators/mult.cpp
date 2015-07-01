#include "function_tree/operators/mult_operator.hpp"

namespace bertini {

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
//            return children_[0];
        }

	
}
