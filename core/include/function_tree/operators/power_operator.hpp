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
// sum_operator.hpp:  Declares the class SumOperator.



#ifndef PowerOperator_h
#define PowerOperator_h

#include <cmath>
#include "function_tree/operators/binary_operator.hpp"
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
		
		virtual dbl FreshEval(dbl) override
		{
			return std::pow( base_->Eval<dbl>(), exponent_->Eval<dbl>());
		}
		
		
		virtual mpfr FreshEval(mpfr) override
		{
			return std::pow( base_->Eval<mpfr>(), exponent_->Eval<mpfr>());
		}
		
	private:
		
		std::shared_ptr<Node> base_;
		std::shared_ptr<Node> exponent_;
	};
	// end of the class PowerOperator
	
	
	
	
	// begin the overload of operators
	
	
	inline std::shared_ptr<bertini::Node> pow(const std::shared_ptr<bertini::Node> & N, const std::shared_ptr<bertini::Node> & p)
	{
		return std::make_shared<bertini::PowerOperator>(N,p);
	}
	
	inline std::shared_ptr<bertini::Node> pow(const std::shared_ptr<bertini::Node> & N, double p)
	{
		return std::make_shared<bertini::PowerOperator>(N,std::make_shared<bertini::Number>(p));
	}
	
	
} // re: namespace bertini


#endif

