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
//  Created by Collins, James B. on 4/30/15.
//
//
// sum_operator.hpp:  Declares the class SumOperator.

#ifndef operator_arithmetic_hpp
#define operator_arithmetic_hpp

#include <vector>
#include <string>

#include "function_tree/node.hpp"

#include "function_tree/operators/operator.hpp"



#include "function_tree/symbols/number.hpp"
#include "function_tree/symbols/variable.hpp"
#include "function_tree/symbols/differential.hpp"

#include <cmath>








namespace bertini {
	
	
	// Node -> NaryOperator -> SumOperator
	// Description: This class represents summation and difference operators.  All children are terms and are stored
	// in a single vector, and a vector of bools is used to determine the sign of each term.  FreshEval method
	// is defined for summation and difference.
	class SumOperator : public virtual NaryOperator
	{
	public:
		virtual ~SumOperator() = default;
		
		
		
		SumOperator(){}
		
		SumOperator(const std::shared_ptr<Node> & left, const std::shared_ptr<Node> & right)
		{
			AddChild(left);
			AddChild(right);
		}
		
		
		SumOperator(const std::shared_ptr<Node> & left, bool add_or_sub_left, const std::shared_ptr<Node> & right, bool add_or_sub_right)
		{
			AddChild(left, add_or_sub_left);
			AddChild(right, add_or_sub_right);
		}
		
		
		SumOperator& operator+=(const std::shared_ptr<Node> & rhs)
		{
			this->AddChild(rhs);
			return *this;
		}
		
		SumOperator& operator-=(const std::shared_ptr<Node> & rhs)
		{
			this->AddChild(rhs,false);
			return *this;
		}
		
		
		
		
		
		
		//Special Behaviour: by default all terms added are positive
		void AddChild(std::shared_ptr<Node> child) override
		{
			NaryOperator::AddChild(std::move(child));
			children_sign_.push_back(true);
		}
		
		
		
		//Special Behaviour: Pass bool to set sign of term: true = add, false = subtract
		void AddChild(std::shared_ptr<Node> child, bool sign) // not an override
		{
			NaryOperator::AddChild(std::move(child));
			children_sign_.push_back(sign);
		}
		
		
		/**
		 Method for printing to output stream
		 */
		void print(std::ostream & target) const override;
		
		
		
		
		/**
		 Return SumOperator whose children are derivatives of children_
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		/**
		 Compute the degree of a node.  For sum functions, the degree is the max among summands.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		int Degree(VariableGroup const& vars) const;

		/**
		 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
	    */
		std::vector<int> MultiDegree(VariableGroup const& vars) const override;
		




		/**
		 Homogenize a sum, with respect to a variable group, and using a homogenizing variable.
		 */
		void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar);
		
		bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override;

		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& vars) const override;
		

	 



	protected:
		/**
		 Specific implementation of FreshEval for add and subtract.
		 If child_sign_ = true, then add, else subtract
		 */
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override;
		
		
		/**
		 Specific implementation of FreshEval for add and subtract.
		 If child_sign_ = true, then add, else subtract
		 */
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override;
		
	private:
		// Stores the sign of the particular term.  There is a one-one
		// correspondence between elements of children_sign_ and children_.  This
		// is enforced by the AddChild method below, redefined in SumOperator.
		
		// TODO(JBC): If we add method to delete child, must also delete children_sign_ entry.
		std::vector<bool> children_sign_;
		
	};
	
	
	
	///////////////////
	//
	//     SUM AND DIFFERENCE ARITHMETIC OPERATORS
	//
	/////////////////////
	
	
	
	///////////////
	//
	//  addition operators
	//
	///////////////
	
	
	inline std::shared_ptr<Node>& operator+=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
	{
		std::shared_ptr<Node> temp = std::make_shared<SumOperator>();
		
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(lhs);
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(rhs);
		
		lhs.swap(temp);
		return lhs;
	}
	
	inline std::shared_ptr<Node>& operator+=(std::shared_ptr<Node> & lhs, double rhs)
	{
		std::shared_ptr<Node> temp = std::make_shared<SumOperator>();
		
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(lhs);
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(std::make_shared<Float>(rhs));
		
		lhs.swap(temp);
		return lhs;
	}
	
	
	
	
	
	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return std::make_shared<SumOperator>(lhs,rhs);
	}
	
	
	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, double rhs)
	{
		return std::make_shared<SumOperator>(lhs,std::make_shared<Float>(rhs));
	}
	
	inline std::shared_ptr<Node> operator+(double lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Float>(lhs), rhs);
	}
	
	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, dbl rhs)
	{
		return std::make_shared<SumOperator>(lhs, std::make_shared<Float>(rhs));
	}
	
	inline std::shared_ptr<Node> operator+(dbl lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Float>(lhs), rhs);
	}

	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, mpfr rhs)
	{
		return std::make_shared<SumOperator>(lhs,std::make_shared<Float>(rhs));
	}
	
	inline std::shared_ptr<Node> operator+(mpfr lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Float>(lhs), rhs);
	}

	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, int rhs)
	{
		return std::make_shared<SumOperator>(lhs,std::make_shared<Integer>(rhs));
	}
	
	inline std::shared_ptr<Node> operator+(int lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Integer>(lhs), rhs);
	}

	
	
	
	
	///////////////
	//
	//  subtraction operators
	//
	///////////////
	
	
	inline std::shared_ptr<Node>& operator-=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
	{
		std::shared_ptr<Node> temp = std::make_shared<SumOperator>();
		
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(lhs);
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(rhs,false);
		
		lhs.swap(temp);
		return lhs;
	}
	
	inline std::shared_ptr<Node>& operator-=(std::shared_ptr<Node> & lhs, double rhs)
	{
		std::shared_ptr<Node> temp = std::make_shared<SumOperator>();
		
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(lhs);
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(std::make_shared<Float>(rhs),false);
		
		lhs.swap(temp);
		return lhs;
	}
	
	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return std::make_shared<SumOperator>(lhs,true,rhs,false);
	}
	
	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, double rhs)
	{
		return std::make_shared<SumOperator>(lhs, true, std::make_shared<Float>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator-(double lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Float>(lhs), true, rhs, false);
	}
	
	
	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, dbl rhs)
	{
		return std::make_shared<SumOperator>(lhs, true, std::make_shared<Float>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator-(dbl lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Float>(lhs), true, rhs, false);
	}

	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, mpfr rhs)
	{
		return std::make_shared<SumOperator>(lhs, true, std::make_shared<Float>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator-(mpfr lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Float>(lhs), true, rhs, false);
	}

	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, int rhs)
	{
		return std::make_shared<SumOperator>(lhs, true, std::make_shared<Integer>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator-(int lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Integer>(lhs), true, rhs, false);
	}

	
	
	
	
	
	
	
	
	
	/**
	 Node -> UnaryOperator -> NegateOperator
	 Description: This class represents the negation operator.  FreshEval method
	 is defined for negation and multiplies the value by -1.
	 */
	class NegateOperator : public virtual UnaryOperator
	{
	public:
		
		NegateOperator(){}
		
		NegateOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		
		
		/**
		 Print to an arbitrary ostream.
		 */
		void print(std::ostream & target) const override;
		
		
		/**
		 Returns negative of derivative of child.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return child_->IsHomogeneous(v);
		}

		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& vars) const override
		{
			return child_->IsHomogeneous(vars);
		}

		virtual ~NegateOperator() = default;
		
	protected:
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override;
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override;
		
	};
	
	
	inline std::shared_ptr<Node> operator-(const std::shared_ptr<Node> & rhs)
	{
		return std::make_shared<NegateOperator>(rhs);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	// Node -> NaryOperator -> MultOperator
	// Description: This class represents multiplication operator.  All children are factors and are stored
	// in a vector.  FreshEval method is defined for multiplication.
	class MultOperator : public virtual NaryOperator
	{
	public:
		//
		
		MultOperator(){}
		
		MultOperator(std::shared_ptr<Node> const& left, std::shared_ptr<Node> const& right)
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
		
		
		
		
		//Special Behaviour: by default all factors are in numerator
		void AddChild(std::shared_ptr<Node> child) override
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
		void print(std::ostream & target) const override;
		
		/**
		 Differentiates using the product rule.  If there is division, consider as ^(-1) and use chain rule.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		
		/**
		 Compute the degree of a node.  For trig functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return nan.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		int Degree(VariableGroup const& vars) const override;

		/**
		 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
	    */
		std::vector<int> MultiDegree(VariableGroup const& vars) const override;
		

		void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override;
		
		bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override;

		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& vars) const override;
	protected:
		
		// Specific implementation of FreshEval for mult and divide.
		//  If child_mult_ = true, then multiply, else divide
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override;
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override;
		
		
		
		
		
		
	private:
		// Stores the mult/div of a factor.  There is a one-one
		// correspondence between elements of children_sign_ and children_.  This
		// is enforced by the AddChild method, redefined in MultOperator.
		
		// TODO(JBC): If we add method to delete child, must also delete children_mult_ entry.
		std::vector<bool> children_mult_or_div_;
		
	};
	
	
	
	
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
		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(std::make_shared<Float>(rhs));
		
		lhs.swap(temp);
		return lhs;
	}
	
	
	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, double rhs)
	{
		return std::make_shared<MultOperator>(lhs,std::make_shared<Float>(rhs));
	}
	
	inline std::shared_ptr<Node> operator*(double lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Float>(lhs), rhs);
	}
	

	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, dbl rhs)
	{
		return std::make_shared<MultOperator>(lhs,std::make_shared<Float>(rhs));
	}
	
	inline std::shared_ptr<Node> operator*(dbl lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Float>(lhs), rhs);
	}

	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, mpfr rhs)
	{
		return std::make_shared<MultOperator>(lhs,std::make_shared<Float>(rhs));
	}
	
	inline std::shared_ptr<Node> operator*(mpfr lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Float>(lhs), rhs);
	}

	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, int rhs)
	{
		return std::make_shared<MultOperator>(lhs,std::make_shared<Integer>(rhs));
	}
	
	inline std::shared_ptr<Node> operator*(int lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Integer>(lhs), rhs);
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
		std::dynamic_pointer_cast<MultOperator>(temp)->AddChild(std::make_shared<Float>(rhs),false);
		
		lhs.swap(temp);
		return lhs;
	}
	
	
	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return lhs/=rhs;
	}
	



	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, double rhs)
	{
		return std::make_shared<MultOperator>(lhs, true, std::make_shared<Float>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator/(double lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Float>(lhs), true, rhs, false);
	}
	

	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, dbl rhs)
	{
		return std::make_shared<MultOperator>(lhs, true, std::make_shared<Float>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator/(dbl lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Float>(lhs), true, rhs, false);
	}

	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, mpfr rhs)
	{
		return std::make_shared<MultOperator>(lhs, true, std::make_shared<Float>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator/(mpfr lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Float>(lhs), true, rhs, false);
	}

	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, int rhs)
	{
		return std::make_shared<MultOperator>(lhs, true, std::make_shared<Integer>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator/(int lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Integer>(lhs), true, rhs, false);
	}
	
	
	/**
	 Operator for power functions with arbitrary expressions in the exponent and base.
	 
	 \see IntegerPowerOperator
	 */
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
		
		
		void Reset() override;
		
		
		
		void print(std::ostream & target) const override;
		
		
		
		/**
		 Differentiates with the power rule.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		
		/**
		 Compute the degree of a node.  For power functions, the degree depends on the degree of the power.  If the exponent is constant, then the degree is actually a number.  If the exponent is non-constant, then the degree is ill-defined.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		int Degree(VariableGroup const& vars) const override;

		/**
		 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
	    */
		std::vector<int> MultiDegree(VariableGroup const& vars) const override;
		


		void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override;
		
		bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override;

		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& vars) const override;

		virtual ~PowerOperator() = default;
		
		/**
		 Change the precision of this variable-precision tree node.
		 
		 \param prec the number of digits to change precision to.
		 */
		virtual void precision(unsigned int prec) override
		{
			auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
			val_pair.first.precision(prec);

			base_->precision(prec);
			exponent_->precision(prec);
		}



	protected:
		
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override;
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override;
		
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
		return std::make_shared<bertini::PowerOperator>(N,std::make_shared<bertini::Float>(p));
	}
	

	inline std::shared_ptr<bertini::Node> pow(const std::shared_ptr<bertini::Node> & N, dbl p)
	{
		return std::make_shared<bertini::PowerOperator>(N,std::make_shared<bertini::Float>(p));
	}

	inline std::shared_ptr<bertini::Node> pow(const std::shared_ptr<bertini::Node> & N, mpfr p)
	{
		return std::make_shared<bertini::PowerOperator>(N,std::make_shared<bertini::Float>(p));
	}
	
	
	
	
	
	
	
	
	
	
	/**
	 Node -> UnaryOperator -> IntegerPowerOperator
	 
	 Description: This class represents the exponentiation operator.  The base is stored in
	 children_, and an extra variable(exponent_) stores the exponent.  FreshEval is
	 defined as the exponention operation.
	 */
	class IntegerPowerOperator : public virtual UnaryOperator
	{
	public:
		
		
		
		
		/**
		 polymorphic method for printing to an arbitrary stream.
		 */
		void print(std::ostream & target) const override;
		
		
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
		std::shared_ptr<Node> Differentiate() override;
		
		
		/**
		 Compute the degree of a node.  For integer power functions, the degree is the product of the degree of the argument, and the power.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		
		bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return child_->IsHomogeneous(v);
		}
		

		/**
		Check for homogeneity, with respect to a variable group.
		*/
		bool IsHomogeneous(VariableGroup const& vars) const override
		{
			return child_->IsHomogeneous(vars);
		}


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
	}; // re: class IntegerPowerOperator
	
	
	
	inline std::shared_ptr<Node> pow(std::shared_ptr<Node> const& base, int power)
	{
		return std::make_shared<bertini::IntegerPowerOperator>(base,power);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 Node -> UnaryOperator -> SqrtOperator
	 Description: This class represents the square root function.  FreshEval method
	 is defined for square root and takes the square root of the child node.
	 */
	class SqrtOperator : public  virtual UnaryOperator
	{
	public:
		
		SqrtOperator(){}
		
		SqrtOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the square root function.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For the square root function, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		virtual ~SqrtOperator() = default;
		
	protected:
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override;
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override;
		
	};
	
	
	
	inline std::shared_ptr<bertini::Node> sqrt(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::SqrtOperator>(N);
	}
	
	
	
	
	// Node -> UnaryOperator -> ExpOperator
	// Description: This class represents the exponential function.  FreshEval method
	// is defined for exponential and takes the exponential of the child node.
	class ExpOperator : public  virtual UnaryOperator
	{
	public:
		
		ExpOperator(){}
		
		ExpOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
	 
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the exponential function.
		 */
		std::shared_ptr<Node> Differentiate() override;
		
		
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		virtual ~ExpOperator() = default;
		
	protected:
		
		// Specific implementation of FreshEval for exponentiate.
		dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override;
		
		mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override;
		
	};
	
	
	inline std::shared_ptr<bertini::Node> exp(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::ExpOperator>(N);
	}
	
	
} // re: namespace bertini






#endif
