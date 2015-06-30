//This file is part of Bertini 2.0.
//
//mult_operator.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mult_operator.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with mult_operator.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// mult_operator.hpp:  Declares the class MultOperator.



#ifndef b2Test_MultOperator_h
#define b2Test_MultOperator_h

#include "function_tree/node.hpp"
#include "function_tree/operators/nary_operator.hpp"
#include "function_tree/operators/sum_operator.hpp"


namespace bertini {

	// Node -> NaryOperator -> MultOperator
	// Description: This class represents multiplication operator.  All children are factors and are stored
	// in a vector.  FreshEval method is defined for multiplication.
	class MultOperator : public virtual NaryOperator
	{
	public:
		//

		MultOperator(){}

		MultOperator(const std::shared_ptr<Node> & left, const std::shared_ptr<Node> & right)
		{
			AddChild(left);
			AddChild(right);
		}


		MultOperator(const std::shared_ptr<Node> & left, bool mult_or_div_left, const std::shared_ptr<Node> & right, bool mult_or_div_right)
		{
			AddChild(left, mult_or_div_left);
			AddChild(right, mult_or_div_right);
		}








		virtual ~MultOperator() = default;

		// See node.hpp for description
		// TODO:(JBC) Implement this for multiplication
		virtual std::string PrintNode() override {return "";}


		//Special Behaviour: by default all factors are in numerator
		virtual void AddChild(std::shared_ptr<Node> child) override
		{
			NaryOperator::AddChild(std::move(child));
			children_mult_or_div_.push_back(true);
		}



		//Special Behaviour: Pass bool to set sign of term: true = mult, false = divide
		void AddChild(std::shared_ptr<Node> child, bool mult) // not an override
		{
			NaryOperator::AddChild(std::move(child));
			children_mult_or_div_.push_back(mult);
		}


		/**
		 overridden method for printing to an output stream
		 */
		virtual void print(std::ostream & target) const override
		{
			target << "(";
			for (auto iter = children_.begin(); iter!= children_.end(); iter++) {
				if (iter==children_.begin())
					if (! *children_mult_or_div_.begin()  )
						target << "1/";  
				(*iter)->print(target);
				if (iter!=(children_.end()-1)){
					if (*(children_mult_or_div_.begin() + (iter-children_.begin())+1)) {
						target << "*";
					}
					else{
						target << "/";
					}

				}
			}
			target << ")";
		}




        /**
         Differentiates using the product rule.  If there is division, consider as ^(-1) and use chain rule.
         */
        virtual std::shared_ptr<Node> Differentiate() override
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




	protected:
		// Specific implementation of FreshEval for mult and divide.
		//  If child_mult_ = true, then multiply, else divide
//		virtual dbl FreshEval(dbl) override
//		{
//			dbl retval{1};
//			for(int ii = 0; ii < children_.size(); ++ii)
//			{
//				if(children_mult_or_div_[ii])
//				{
//					retval *= children_[ii]->Eval<dbl>();
//				}
//				else
//				{
//					retval /= children_[ii]->Eval<dbl>();
//				}
//			}
//
//			return retval;
//		}
//
//
//
//
//		virtual mpfr FreshEval(mpfr) override
//		{
//			mpfr retval{1};
//			for(int ii = 0; ii < children_.size(); ++ii)
//			{
//				if(children_mult_or_div_[ii])
//				{
//					retval *= children_[ii]->Eval<mpfr>();
//				}
//				else
//				{
//					retval /= children_[ii]->Eval<mpfr>();
//				}
//			}
//
//			return retval;
//		}


        virtual dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
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




        virtual mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
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






	private:
		// Stores the mult/div of a factor.  There is a one-one
		// correspondence between elements of children_sign_ and children_.  This
		// is enforced by the AddChild method, redefined in MultOperator.

		// TODO(JBC): If we add method to delete child, must also delete children_mult_ entry.
		std::vector<bool> children_mult_or_div_;

	};


} // re: namespace bertini


namespace {


	using Node = bertini::Node;
	using MultOperator = bertini::MultOperator;

	/*
	multiplication operators
	*/
	inline std::shared_ptr<Node> operator*=(std::shared_ptr<MultOperator> & lhs, const std::shared_ptr<Node> & rhs)
	{
		lhs->AddChild(rhs);
		return lhs;
	}


	inline std::shared_ptr<Node>& operator*=(std::shared_ptr<Node> & lhs, double rhs)
	{
		std::shared_ptr<Node> temp = std::make_shared<MultOperator>();

		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(lhs);
		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(std::make_shared<Number>(rhs));

		lhs.swap(temp);
		return lhs;
	}

	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, double rhs)
	{
		return std::make_shared<MultOperator>(lhs,std::make_shared<Number>(rhs));
	}

	inline std::shared_ptr<Node> operator*(double lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Number>(lhs), rhs);
	}

	inline std::shared_ptr<Node> operator*(int lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Number>(lhs), rhs);
	}

	inline std::shared_ptr<Node>& operator*=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
	{

		std::shared_ptr<Node> temp = std::make_shared<MultOperator>();

		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(lhs);
		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(rhs);

		lhs.swap(temp);
		return lhs;
	}


	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return lhs*=rhs;
	}






	/*
	division operators
	*/


	inline std::shared_ptr<Node>& operator/=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
	{

		std::shared_ptr<Node> temp = std::make_shared<MultOperator>();

		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(lhs);
		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(rhs,false);

		lhs.swap(temp);
		return lhs;
	}


	inline std::shared_ptr<Node> operator/=(std::shared_ptr<MultOperator> & lhs, const std::shared_ptr<Node> & rhs)
	{
		lhs->AddChild(rhs,false);
		return lhs;
	}


	inline std::shared_ptr<Node>& operator/=(std::shared_ptr<Node> & lhs, double rhs)
	{
		std::shared_ptr<Node> temp = std::make_shared<MultOperator>();

		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(lhs);
		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(std::make_shared<Number>(rhs),false);

		lhs.swap(temp);
		return lhs;
	}


	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return lhs/=rhs;
	}

	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, double rhs)
	{
		return std::make_shared<MultOperator>(lhs, true, std::make_shared<Number>(rhs), false);
	}

	inline std::shared_ptr<Node> operator/(double lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Number>(lhs), true, rhs, false);
	}



}



#endif
