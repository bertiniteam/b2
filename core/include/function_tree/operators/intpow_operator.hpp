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
// exp_operator.hpp:  Declares the class  IntegerPowerOperator.



#ifndef b2Test_IntegerPowerOperator_h
#define b2Test_IntegerPowerOperator_h

#include "function_tree/Node.hpp"
#include "function_tree/operators/unary_operator.hpp"




namespace bertini {

	




	using std::pow;

	/**
	 Node -> UnaryOperator -> IntegerPowerOperator

	 Description: This class represents the exponentiation operator.  The base is stored in
	 children_, and an extra variable(exponent_) stores the exponent.  FreshEval is
	 defined as the exponention operation.
	 */
	class IntegerPowerOperator : public virtual UnaryOperator
	{
	public:
		// These do nothing for a constant
		std::string PrintNode() override
		{
			return child_->PrintNode() + "^" + std::to_string(exponent());
		}



		/**
		 Virtual polymorphic method for printing to an arbitrary stream.
		 */
		virtual void print(std::ostream & target) const override;


		/**
		 Get the integet exponent of an ExpOperator
		 */
		void set_exponent(int exp)
		{
			exponent_ = exp;
		}


		/**
		 Get the exponent of an ExpOperator
		 */
		int exponent() const
		{
			return exponent_;
		}


        /**
         Differentiates a number.
         */
        virtual std::shared_ptr<Node> Differentiate() override;





		 /**
		Compute the degree of a node.  For integer power functions, the degree is the product of the degree of the argument, and the power.
        */
		virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;

		


		virtual ~IntegerPowerOperator() = default;


		/**
		 Constructor, passing in the Node you want as the base, and the integer you want for the power.
		 */
		IntegerPowerOperator(const std::shared_ptr<Node> & N, int p = 1) : exponent_(p), UnaryOperator(N)
		{}


		IntegerPowerOperator(){}

	protected:


        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            return pow(child_->Eval<dbl>(diff_variable), exponent_);
        }

        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            return pow(child_->Eval<mpfr>(diff_variable),exponent_);
        }






	private:

		int exponent_ = 1; ///< Exponent for the exponenetial operator
	};



	inline std::shared_ptr<Node> pow(std::shared_ptr<Node> const& base, int power)
	{
		return std::make_shared<bertini::IntegerPowerOperator>(base,power);
	}

}// re: namespace bertini




#endif
