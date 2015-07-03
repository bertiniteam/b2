//This file is part of Bertini 2.0.
//
//trig.hpp is free software: you can redistribute it and/or modify
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
// trig.hpp:  Declares the trigonometric operator classes.


#ifndef trig_operator_hpp
#define trig_operator_hpp

#include "function_tree/operators/operator.hpp"
#include "function_tree/operators/arithmetic.hpp"



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
		
		
		
		
		void print(std::ostream & target) const override
		{
			target << "sin(";
			child_->print(target);
			target << ")";
		}
		
		
		
		/**
		 Differentiates the sine function.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			if (child_->Degree(v)==0)
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
		
		
		
		virtual ~SinOperator() = default;
		
	protected:
		
		
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
	

	// Node -> UnaryOperator -> ArcSinOperator
	// Description: This class represents the inverse sine function.  FreshEval method
	// is defined for arcsine and takes the sine of the child node.
	class ArcSinOperator : public  virtual UnaryOperator
	{
	public:
		
		ArcSinOperator(){}
		
		ArcSinOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		
		
		
		void print(std::ostream & target) const override
		{
			target << "sin(";
			child_->print(target);
			target << ")";
		}
		
		
		
		/**
		 Differentiates the sine function.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			if (child_->Degree(v)==0)
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
		
		
		
		virtual ~ArcSinOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
		{
			return asin(child_->Eval<dbl>(diff_variable));
		}
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
		{
			return asin(child_->Eval<mpfr>(diff_variable));
		}
		
	};


	// begin the overload of operators
	inline std::shared_ptr<bertini::Node> sin(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::SinOperator>(N);
	}
	
	inline std::shared_ptr<bertini::Node> asin(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::ArcSinOperator>(N);
	}
	
	
	
	
	// Node -> UnaryOperator -> CosOperator
	// Description: This class represents the cosine function.  FreshEval method
	// is defined for cosine and takes the cosine of the child node.
	class CosOperator : public  virtual UnaryOperator
	{
	public:
		
		CosOperator(){}
		
		CosOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		
		
		
		void print(std::ostream & target) const override;
		
		
		
		
		/**
		 Differentiates the cosine function.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For trig functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		virtual ~CosOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
		{
			return cos(child_->Eval<dbl>(diff_variable));
		}
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
		{
			return cos(child_->Eval<mpfr>(diff_variable));
		}
		
	};
	
	
	// Node -> UnaryOperator -> ArcCosOperator
	// Description: This class represents the inverse cosine function.  FreshEval method
	// is defined for arccosine and takes the arccosine of the child node.
	class ArcCosOperator : public  virtual UnaryOperator
	{
	public:
		
		ArcCosOperator(){}
		
		ArcCosOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		
		
		
		void print(std::ostream & target) const override;
		
		
		
		
		/**
		 Differentiates the cosine function.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For trig functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		virtual ~ArcCosOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
		{
			return acos(child_->Eval<dbl>(diff_variable));
		}
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
		{
			return acos(child_->Eval<mpfr>(diff_variable));
		}
		
	};
	
	
	inline std::shared_ptr<bertini::Node> cos(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::CosOperator>(N);
	}



	inline std::shared_ptr<bertini::Node> acos(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::ArcCosOperator>(N);
	}
	
	


	
	
	
	// Node -> UnaryOperator -> TanOperator
	// Description: This class represents the tangent function.  FreshEval method
	// is defined for tangent and takes the tangent of the child node.
	class TanOperator : public  virtual UnaryOperator
	{
	public:
		
		TanOperator(){}
		
		TanOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		
		
		
		void print(std::ostream & target) const override
		{
			target << "tan(";
			child_->print(target);
			target << ")";
		}
		
		
		/**
		 Differentiates the tangent function.
		 */
		std::shared_ptr<Node> Differentiate() override
		{
			auto ret_mult = std::make_shared<MultOperator>();
			ret_mult->AddChild(child_->Differentiate());
			ret_mult->AddChild(std::make_shared<CosOperator>(child_),false);
			ret_mult->AddChild(std::make_shared<CosOperator>(child_),false);
			return ret_mult;
		}
		
		
		
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			if (child_->Degree(v)==0)
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
		
		
		virtual ~TanOperator() = default;
		
	protected:
		
		
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
	
	
	// Node -> UnaryOperator -> ArcTanOperator
	// Description: This class represents the inverse tangent function.  FreshEval method
	// is defined for arctangent and takes the arc tangent of the child node.
	class ArcTanOperator : public  virtual UnaryOperator
	{
	public:
		
		ArcTanOperator(){}
		
		ArcTanOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		
		
		
		void print(std::ostream & target) const override
		{
			target << "atan(";
			child_->print(target);
			target << ")";
		}
		
		
		/**
		 Differentiates the tangent function.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			if (child_->Degree(v)==0)
			{
				return 0;
			}
			else
			{
				return -1;
			}
		}
		
		
		virtual ~ArcTanOperator() = default;
		
	protected:
		
		
		// Specific implementation of FreshEval for arctangent.
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
		{
			return atan(child_->Eval<dbl>(diff_variable));
		}
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
		{
			return atan(child_->Eval<mpfr>(diff_variable));
		}
		
	};
	
	
	// begin the overload of operators
	
	
	inline std::shared_ptr<bertini::Node> tan(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::TanOperator>(N);
	}


	
	inline std::shared_ptr<bertini::Node> atan(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::ArcTanOperator>(N);
	}
	
	
	
} // re: namespace bertini


#endif
