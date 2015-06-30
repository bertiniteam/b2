//This file is part of Bertini 2.0.
//
//sum_operator.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//sum_operator.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with sum_operator.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Daniel A Brake
//
//
// power_operator.hpp:  Declares the class PowerOperator.



#ifndef PowerOperator_h
#define PowerOperator_h

#include <cmath>
#include "function_tree/operators/binary_operator.hpp"
#include "function_tree/operators/mult_operator.hpp"

#include "function_tree/symbols/number.hpp"


namespace bertini {

	class PowerOperator : public virtual BinaryOperator
	{

	public:

		PowerOperator(){}

		PowerOperator(const std::shared_ptr<Node> & new_base, const std::shared_ptr<Node> & new_exponent) : base_(new_base), exponent_(new_exponent)
		{
		}



		void SetBase(std::shared_ptr<Node> new_base)
		{
			base_ = new_base;
		}

		void SetExponent(std::shared_ptr<Node> new_exponent)
		{
			exponent_ = new_exponent;
		}







		virtual void Reset()
		{
			Node::ResetStoredValues();
			base_->Reset();
			exponent_->Reset();
		}



		virtual void print(std::ostream & target) const override
		{
			target << "(" << *base_ << ")^(" << *exponent_ << ")";
		}



        /**
         Differentiates with the power rule.
         */
        virtual std::shared_ptr<Node> Differentiate() override
        {
            auto ret_mult = std::make_shared<MultOperator>();
            auto exp_minus_one = std::make_shared<SumOperator>(exponent_, true, std::make_shared<Number>("1.0"),false);
            ret_mult->AddChild(base_->Differentiate());
            ret_mult->AddChild(exponent_);
            ret_mult->AddChild(std::make_shared<PowerOperator>(base_, exp_minus_one));
            return ret_mult;
        }


        


		 /**
		Compute the degree of a node.  For power functions, the degree depends on the degree of the power.  If the exponent is constant, then the degree is actually a number.  If the exponent is non-constant, then the degree is ill-defined.
        */
		virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
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



		virtual ~PowerOperator() = default;






		////////////// TESTING /////////////////
		virtual void PrintTree() override
		{
			for(int ii = 0; ii < tabcount; ++ii)
			{
				std::cout << "\t";
			}
			std::cout << tabcount+1 << "." <<  boost::typeindex::type_id_runtime(*this).pretty_name() << " = " << this->Eval<dbl>()<< std::endl;
			tabcount++;
			base_->PrintTree();
			exponent_->PrintTree();
			tabcount--;
		}


		// Print the data for this node, as well as all it's children
		//TODO (JBC): Implement this method
		virtual std::string PrintNode() override {return "";}
		////////////// TESTING /////////////////




	protected:

//		virtual dbl FreshEval(dbl) override
//		{
//			return std::pow( base_->Eval<dbl>(), exponent_->Eval<dbl>());
//		}
//
//
//		virtual mpfr FreshEval(mpfr) override
//		{
//			return pow( base_->Eval<mpfr>(), exponent_->Eval<mpfr>());
//		}


        virtual dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            return std::pow( base_->Eval<dbl>(diff_variable), exponent_->Eval<dbl>());
        }


        virtual mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            return pow( base_->Eval<mpfr>(diff_variable), exponent_->Eval<mpfr>());
        }

	private:

		std::shared_ptr<Node> base_;
		std::shared_ptr<Node> exponent_;
	};
	// end of the class PowerOperator






} // re: namespace bertini



// begin the overload of operators


	inline std::shared_ptr<bertini::Node> pow(const std::shared_ptr<bertini::Node> & N, const std::shared_ptr<bertini::Node> & p)
	{
		return std::make_shared<bertini::PowerOperator>(N,p);
	}

	inline std::shared_ptr<bertini::Node> pow(const std::shared_ptr<bertini::Node> & N, double p)
	{
		return std::make_shared<bertini::PowerOperator>(N,std::make_shared<bertini::Number>(p));
	}





#endif

