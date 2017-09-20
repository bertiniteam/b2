//This file is part of Bertini 2.
//
//arithmetic.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//arithmetic.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with arithmetic.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
// Dani Brake
// University of Notre Dame
//
//  Created by Collins, James B. on 4/30/15.
//
//
// arithmetic.hpp:  Declares the arithmetic nodes for bertini2.

/**
\file arithmetic.hpp

\brief Provides the artithmetic Node types, such as Sum and Power

*/

#ifndef BERTINI_OPERATOR_ARITHMETIC_HPP
#define BERTINI_OPERATOR_ARITHMETIC_HPP

#include <vector>
#include <string>

#include "bertini2/function_tree/node.hpp"

#include "bertini2/function_tree/operators/operator.hpp"



#include "bertini2/function_tree/symbols/number.hpp"
#include "bertini2/function_tree/symbols/special_number.hpp"
#include "bertini2/function_tree/symbols/variable.hpp"
#include "bertini2/function_tree/symbols/differential.hpp"

#include "bertini2/function_tree/forward_declares.hpp"

#include <cmath>








namespace bertini {
namespace node{	
	/**
	\brief Represents summation and difference Operator.

	This class represents summation and difference operators.  All children are terms and are stored
	in a single vector, and a vector of bools is used to determine the sign of each term.  FreshEval method
	is defined for summation and difference.
	*/
	class SumOperator : public virtual NaryOperator
	{
	public:
		virtual ~SumOperator() = default;
		
		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;
		unsigned ReduceDepth() override;
		unsigned ReduceSubSums();
		unsigned ReduceSubMults();

		SumOperator(const std::shared_ptr<Node> & s, bool add_or_sub)
		{
			AddChild(s, add_or_sub);
		}
		
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
		
		
		
		
		
		/**
		\note: Special Behaviour: by default all terms added are positive
		*/
		void AddChild(std::shared_ptr<Node> child) override
		{
			NaryOperator::AddChild(std::move(child));
			children_sign_.push_back(true);
		}
		
		
		/**
		\note Special Behaviour: Pass bool to set sign of term: true = add, false = subtract
		*/
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
		 Return SumOperator whose children are derivatives of the children, omitted as possible
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;

		/**
		 Compute the degree of a node.  For sum functions, the degree is the max among summands.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		int Degree(VariableGroup const& vars) const override;

		/**
		 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
		*/
		std::vector<int> MultiDegree(VariableGroup const& vars) const override;
		




		/**
		 Homogenize a sum, with respect to a variable group, and using a homogenizing variable.
		 */
		void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override;
		
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
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;

		/**
		 Specific implementation of FreshEval in place for add and subtract.
		 If child_sign_ = true, then add, else subtract
		 */
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		
		/**
		 Specific implementation of FreshEval for add and subtract.
		 If child_sign_ = true, then add, else subtract
		 */
		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;

		/**
		 Specific implementation of FreshEval for add and subtract.
		 If child_sign_ = true, then add, else subtract
		 */
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

	private:
		// Stores the sign of the particular term.  There is a one-one
		// correspondence between elements of children_sign_ and children_.  This
		// is enforced by the AddChild method below, redefined in SumOperator.
		
		// TODO(JBC): If we add method to delete child, must also delete children_sign_ entry.
		std::vector<bool> children_sign_;
		
	private:

		SumOperator() = default;

		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<NaryOperator>(*this);
			ar & children_sign_;
		}


		void PrecisionChangeSpecific(unsigned prec) const override
		{
			temp_mp_.precision(prec);
		}

		mutable mpfr temp_mp_;
		mutable dbl temp_d_;

		friend MultOperator;
	};
	
	
	

	
	
	
	
	
	
	
	
	
	/**
	\brief The negation Operator.

	 This class represents the negation Operator.  FreshEval method
	 is defined for negation and multiplies the value by -1.
	 */
	class NegateOperator : public virtual UnaryOperator
	{
	public:
		
		NegateOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		
		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;

		/**
		 Print to an arbitrary ostream.
		 */
		void print(std::ostream & target) const override;
		
		
		/**
		 Returns negative of derivative of child.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;

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
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


	private:

		NegateOperator() = default;

		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<UnaryOperator>(*this);
		}
	};
	
	
	inline std::shared_ptr<Node> operator-(const std::shared_ptr<Node> & rhs)
	{
		return std::make_shared<NegateOperator>(rhs);
	}
	
	
	
	
	
	
	
	
	
	
	
	/**
	\brief Multiplication and division Operator.

	This class represents the Operator for multiplication and division.  All children are factors and are stored
	in a vector.  FreshEval method is defined for multiplication.
	*/
	class MultOperator : public virtual NaryOperator
	{
	public:

		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;
		unsigned ReduceDepth() override;
		unsigned ReduceSubSums();
		unsigned ReduceSubMults();

		/**
		 single-node instantiation.  
		
		 if evaluated, would return simply the value of s.
		 */
		MultOperator(std::shared_ptr<Node> const& s)
		{
			AddChild(s);
		}
		
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
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
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
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		
		
		
		
		
	private:
		
		MultOperator() = default;

		// Stores the mult/div of a factor.  There is a one-one
		// correspondence between elements of children_sign_ and children_.  This
		// is enforced by the AddChild method, redefined in MultOperator.
		
		// TODO(JBC): If we add method to delete child, must also delete children_mult_ entry.
		std::vector<bool> children_mult_or_div_;
		

	private:

		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<NaryOperator>(*this);
			ar & children_mult_or_div_;
		}

		void PrecisionChangeSpecific(unsigned prec) const override
		{
			temp_mp_.precision(prec);
		}

		mutable mpfr temp_mp_;
		mutable dbl temp_d_;

		friend SumOperator;
	};
	
	
	
	
	
	/**
	\brief Operator for power functions with arbitrary expressions in the exponent and base.

	Operator for power functions with arbitrary expressions in the exponent and base.
	 
	 
	 \see IntegerPowerOperator
	 */
	class PowerOperator : public virtual Operator
	{
		
	public:
		
		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;

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
		
		
		void Reset() const override;
		
		
		
		void print(std::ostream & target) const override;
		
		
		
		/**
		 Differentiates with the power rule.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
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
		virtual void precision(unsigned int prec) const override
		{
			auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
			val_pair.first.precision(prec);

			base_->precision(prec);
			exponent_->precision(prec);
		}



	protected:
		
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_d(dbl& evaulation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_mp(mpfr& evaulation_value, std::shared_ptr<Variable> const& diff_variable) const override;

	private:
				
		PowerOperator() = default;
		
		std::shared_ptr<Node> base_;
		std::shared_ptr<Node> exponent_;
		


	private:

		friend class boost::serialization::access;
		

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Operator>(*this);
			ar & base_;
			ar & exponent_;
		}
	};
	// end of the class PowerOperator
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	\brief This class represents the exponentiation Operator.


	 This class represents the exponentiation operator.  The base is stored in
	 children_, and an extra variable(exponent_) stores the exponent.  FreshEval is
	 defined as the exponention operation.
	 */
	class IntegerPowerOperator : public virtual UnaryOperator
	{
	public:
		
		
		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;
		
		/**
		 polymorphic method for printing to an arbitrary stream.
		 */
		void print(std::ostream & target) const override;
		
		
		/**
		 Set the integer exponent of an integer power operator
		 */
		void set_exponent(int exp)
		{
			exponent_ = exp;
		}
		
		
		/**
		 Get the exponent
		 */
		int exponent() const
		{
			return exponent_;
		}
		
		
		/**
		 \brief Differentiate
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
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
		IntegerPowerOperator(const std::shared_ptr<Node> & N, int p) : exponent_(p), UnaryOperator(N)
		{}
		
		
		
	protected:
		
		
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override
		{
			return pow(child_->Eval<dbl>(diff_variable), exponent_);
		}

		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override
		{
			child_->EvalInPlace<dbl>(evaluation_value, diff_variable);
			evaluation_value = pow(evaluation_value, exponent_);
		}

		
		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override
		{
			return pow(child_->Eval<mpfr>(diff_variable),exponent_);
		}

		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override
		{
			child_->EvalInPlace<mpfr>(evaluation_value, diff_variable);
			evaluation_value = pow(evaluation_value, exponent_);
		}

	private:
		
		IntegerPowerOperator() = default;


		int exponent_; ///< Exponent for the exponenetial operator

		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<UnaryOperator>(*this);
			ar & exponent_;
		}
	}; // re: class IntegerPowerOperator
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	\brief Represents the square root Operator


	 This class represents the square root function.  FreshEval method
	 is defined for square root and takes the square root of the child node.
	 */
	class SqrtOperator : public  virtual UnaryOperator
	{
	public:
		
		SqrtOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		
		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the square root function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For the square root function, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		virtual ~SqrtOperator() = default;
		
	protected:
		
		// Specific implementation of FreshEval for negate.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;


	private:

		SqrtOperator() = default;

		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<UnaryOperator>(*this);
		}
	};
	
	
	
	inline std::shared_ptr<Node> sqrt(const std::shared_ptr<Node> & N)
	{
		return std::make_shared<SqrtOperator>(N);
	}
	
	
	
	
	/**
	\brief represents the exponential function

	This class represents the exponential function.  FreshEval method
	is defined for exponential and takes the exponential of the child node.
	*/
	class ExpOperator : public  virtual UnaryOperator
	{
	public:

		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;

		ExpOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
	 
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the exponential function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		virtual ~ExpOperator() = default;
		
	protected:
		
		// Specific implementation of FreshEval for exponentiate.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;

	private:
		ExpOperator() = default;
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<UnaryOperator>(*this);
		}
	};
	
	
	/**
	\brief represents the natural logarithm function

	This class represents the natural logarithm function.
	*/
	class LogOperator : public  virtual UnaryOperator
	{
	public:
		
		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;


		LogOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
	 
		
		
		
		void print(std::ostream & target) const override;
		
		
		/**
		 Differentiates the exponential function.
		 */
		std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override;
		
		/**
		 Compute the degree with respect to a single variable.
		 
		 For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		 */
		int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
		

		virtual ~LogOperator() = default;
		
	protected:
		
		// Specific implementation of FreshEval for exponentiate.
		dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
		mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override;
		void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override;
		
	private:
		LogOperator() = default;
		friend class boost::serialization::access;
		
		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<UnaryOperator>(*this);
		}
	};

	


	// begin the overload of operators

	inline std::shared_ptr<Node> exp(const std::shared_ptr<Node> & N)
	{
		return std::make_shared<ExpOperator>(N);
	}
	
	inline std::shared_ptr<Node> log(const std::shared_ptr<Node> & N)
	{
		return std::make_shared<LogOperator>(N);
	}
	
	inline std::shared_ptr<Node> pow(const std::shared_ptr<Node> & N, const std::shared_ptr<Node> & p)
	{
		return std::make_shared<PowerOperator>(N,p);
	}

	inline std::shared_ptr<Node> pow(std::shared_ptr<Node> const& base, int power)
	{
		return std::make_shared<IntegerPowerOperator>(base,power);
	}

	std::shared_ptr<Node> pow(const std::shared_ptr<Node> & N, double p) = delete;
	
	std::shared_ptr<Node> pow(const std::shared_ptr<Node> & N, dbl p) = delete;

	inline std::shared_ptr<Node> pow(const std::shared_ptr<Node> & N, mpfr_float p)
	{
		return std::make_shared<PowerOperator>(N,std::make_shared<Float>(p));
	}

	inline std::shared_ptr<Node> pow(const std::shared_ptr<Node> & N, mpfr p)
	{
		return std::make_shared<PowerOperator>(N,std::make_shared<Float>(p));
	}

	inline std::shared_ptr<Node> pow(const std::shared_ptr<Node> & N, mpq_rational const& p)
	{
		return std::make_shared<PowerOperator>(N,std::make_shared<Rational>(p,0));
	}






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
		std::shared_ptr<Node> temp = std::make_shared<SumOperator>(lhs,rhs);		
		lhs.swap(temp);
		return lhs;
	}
		
	
	
	
	
	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return std::make_shared<SumOperator>(lhs,rhs);
	}
	
	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, mpfr_float const& rhs)
	{
		return std::make_shared<SumOperator>(lhs,std::make_shared<Float>(rhs));
	}

	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, mpfr const& rhs)
	{
		return std::make_shared<SumOperator>(lhs,std::make_shared<Float>(rhs));
	}
	
	inline std::shared_ptr<Node> operator+(mpfr_float const& lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Float>(lhs), rhs);
	}

	inline std::shared_ptr<Node> operator+(mpfr const& lhs,  std::shared_ptr<Node> rhs)
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

	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, mpz_int const& rhs)
	{
		return std::make_shared<SumOperator>(lhs,std::make_shared<Integer>(rhs));
	}
	
	inline std::shared_ptr<Node> operator+(mpz_int const& lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Integer>(lhs), rhs);
	}

	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, mpq_rational const& rhs)
	{
		return std::make_shared<SumOperator>(lhs,std::make_shared<Rational>(rhs,0));
	}
	
	inline std::shared_ptr<Node> operator+(mpq_rational const& lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Rational>(lhs,0), rhs);
	}
	
	
	
	///////////////
	//
	//  subtraction operators
	//
	///////////////
	
	
	inline std::shared_ptr<Node>& operator-=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
	{
		std::shared_ptr<Node> temp = std::make_shared<SumOperator>(lhs,true,rhs,false);
		lhs.swap(temp);
		return lhs;
	}
		
	
	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return std::make_shared<SumOperator>(lhs,true,rhs,false);
	}
	
	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, mpfr_float rhs)
	{
		return std::make_shared<SumOperator>(lhs, true, std::make_shared<Float>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator-(mpfr_float lhs,  std::shared_ptr<Node> rhs)
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

	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, mpz_int const& rhs)
	{
		return std::make_shared<SumOperator>(lhs, true, std::make_shared<Integer>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator-(mpz_int const& lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Integer>(lhs), true, rhs, false);
	}

	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, mpq_rational const& rhs)
	{
		return std::make_shared<SumOperator>(lhs, true, std::make_shared<Rational>(rhs,0), false);
	}
	
	inline std::shared_ptr<Node> operator-(mpq_rational const& lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Rational>(lhs,0), true, rhs, false);
	}


	
	/*
	 multiplication operators
	 */
	inline std::shared_ptr<Node> operator*=(std::shared_ptr<MultOperator> & lhs, const std::shared_ptr<Node> & rhs)
	{
		lhs->AddChild(rhs);
		return lhs;
	}
	
	
	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, mpfr_float rhs)
	{
		return std::make_shared<MultOperator>(lhs,std::make_shared<Float>(rhs));
	}

	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, mpfr rhs)
	{
		return std::make_shared<MultOperator>(lhs,std::make_shared<Float>(rhs));
	}
	
	inline std::shared_ptr<Node> operator*(mpfr_float lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Float>(lhs), rhs);
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
	
	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, mpz_int const& rhs)
	{
		return std::make_shared<MultOperator>(lhs,std::make_shared<Integer>(rhs));
	}
	
	inline std::shared_ptr<Node> operator*(mpz_int const& lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Integer>(lhs), rhs);
	}
	
	inline std::shared_ptr<Node> operator*(std::shared_ptr<Node> lhs, mpq_rational const& rhs)
	{
		return std::make_shared<MultOperator>(lhs,std::make_shared<Rational>(rhs,0));
	}
	
	inline std::shared_ptr<Node> operator*(mpq_rational const& lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Rational>(lhs,0), rhs);
	}


	// this function provides an optimization for combining two power operators with the same base.
	inline std::shared_ptr<Node>& operator*=(std::shared_ptr<Node> & lhs, const std::shared_ptr<Node> & rhs)
	{
		
		// if the two nodes are integer power operators, and if they point the same place, then add the powers.

		if (std::dynamic_pointer_cast<IntegerPowerOperator>(lhs) && std::dynamic_pointer_cast<IntegerPowerOperator>(rhs))
		{

			auto lhs_as_intpow = std::dynamic_pointer_cast<IntegerPowerOperator>(lhs);
			auto rhs_as_intpow = std::dynamic_pointer_cast<IntegerPowerOperator>(rhs);
			if (lhs_as_intpow->first_child()==rhs_as_intpow->first_child())
			{
				if (lhs_as_intpow->exponent()>=0 && rhs_as_intpow->exponent()>=0)
				{
					std::shared_ptr<Node> temp = pow(lhs_as_intpow->first_child(),lhs_as_intpow->exponent() + rhs_as_intpow->exponent());
					lhs.swap(temp);
					return lhs;
				}
			}
		}

		std::shared_ptr<Node> temp = std::make_shared<MultOperator>(lhs,rhs);
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

		// if (std::dynamic_pointer_cast<IntegerPowerOperator>(lhs) && std::dynamic_pointer_cast<IntegerPowerOperator>(rhs))
		// {

		// 	auto lhs_as_intpow = std::dynamic_pointer_cast<IntegerPowerOperator>(lhs);
		// 	auto rhs_as_intpow = std::dynamic_pointer_cast<IntegerPowerOperator>(rhs);
		// 	if (lhs_as_intpow->first_child()==rhs_as_intpow->first_child())
		// 	{
		// 		std::shared_ptr<Node> temp = pow(lhs_as_intpow->first_child(),lhs_as_intpow->exponent() - rhs_as_intpow->exponent());
		// 		lhs.swap(temp);
		// 		return lhs;
		// 	}
		// }


		std::shared_ptr<Node> temp = std::make_shared<MultOperator>(lhs,true,rhs,false);
		lhs.swap(temp);
		return lhs;
	}
	
	inline std::shared_ptr<Node> operator/=(std::shared_ptr<MultOperator> & lhs, const std::shared_ptr<Node> & rhs)
	{
		lhs->AddChild(rhs,false);
		return lhs;
	}
	
	
	
	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return lhs/=rhs;
	}
	
	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, mpfr_float rhs)
	{
		return std::make_shared<MultOperator>(lhs, true, std::make_shared<Float>(rhs), false);
	}

	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, mpfr rhs)
	{
		return std::make_shared<MultOperator>(lhs, true, std::make_shared<Float>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator/(mpfr_float lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Float>(lhs), true, rhs, false);
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

	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, mpz_int const& rhs)
	{
		return std::make_shared<MultOperator>(lhs, true, std::make_shared<Integer>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator/(mpz_int const& lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Integer>(lhs), true, rhs, false);
	}

	inline std::shared_ptr<Node> operator/(std::shared_ptr<Node> lhs, mpq_rational const& rhs)
	{
		return std::make_shared<MultOperator>(lhs, true, std::make_shared<Rational>(rhs,0), false);
	}
	
	inline std::shared_ptr<Node> operator/(mpq_rational const& lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<MultOperator>(std::make_shared<Rational>(lhs,0), true, rhs, false);
	}



} // re: namespace node	
} // re: namespace bertini






#endif
