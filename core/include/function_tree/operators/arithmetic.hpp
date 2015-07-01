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
		virtual void AddChild(std::shared_ptr<Node> child) override
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
		virtual void print(std::ostream & target) const override;
        
        
        
        
        /**
         Return SumOperator whose children are derivatives of children_
         */
        virtual std::shared_ptr<Node> Differentiate() override;


		/**
		Compute the degree of a node.  For sum functions, the degree is the max among summands.
        */
		virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;


		/**
		Homogenize a sum, with respect to a variable group, and using a homogenizing variable.
		*/
		void Homogenize(std::vector< std::shared_ptr< Variable > > const& vars, std::shared_ptr<Variable> const& homvar);

		
	protected:
		/**
		Specific implementation of FreshEval for add and subtract.
		 If child_sign_ = true, then add, else subtract
		 */
		virtual dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override;
        
        
        /**
		Specific implementation of FreshEval for add and subtract.
		 If child_sign_ = true, then add, else subtract
		 */
        virtual mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override;
		
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
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(std::make_shared<Number>(rhs));
		
		lhs.swap(temp);
		return lhs;
	}
	
	
	
	
	
	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return std::make_shared<SumOperator>(lhs,rhs);
	}
	
	
	inline std::shared_ptr<Node> operator+(std::shared_ptr<Node> lhs, double rhs)
	{
		return std::make_shared<SumOperator>(lhs,std::make_shared<Number>(rhs));
	}
	
	inline std::shared_ptr<Node> operator+(double lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Number>(lhs), rhs);
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
		std::dynamic_pointer_cast<SumOperator>(temp)->AddChild(std::make_shared<Number>(rhs),false);
		
		lhs.swap(temp);
		return lhs;
	}
	
	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, const std::shared_ptr<Node> & rhs)
	{
		return std::make_shared<SumOperator>(lhs,true,rhs,false);
	}
	
	inline std::shared_ptr<Node> operator-(std::shared_ptr<Node> lhs, double rhs)
	{
		return std::make_shared<SumOperator>(lhs, true, std::make_shared<Number>(rhs), false);
	}
	
	inline std::shared_ptr<Node> operator-(double lhs,  std::shared_ptr<Node> rhs)
	{
		return std::make_shared<SumOperator>(std::make_shared<Number>(lhs), true, rhs, false);
	}
	


	// Node -> UnaryOperator -> NegateOperator
	// Description: This class represents the negation operator.  FreshEval method
	// is defined for negation and multiplies the value by -1.
	class NegateOperator : public  virtual UnaryOperator
	{
	public:
		
		NegateOperator(){}
		
		NegateOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		

		
		
		
		virtual void print(std::ostream & target) const override
		{
			target << "-(";
			child_->print(target);
			target << ")";
		}
        
        
        
        
        /**
         Returns negative of derivative of child.
         */
        virtual std::shared_ptr<Node> Differentiate() override
        {
            return std::make_shared<NegateOperator>(child_->Differentiate());
        }

		 /**
		Compute the degree of a node.  For trig functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return nan.
        */
		virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			return child_->Degree(v);
		}



		virtual ~NegateOperator() = default;
		
	protected:
        
        // Specific implementation of FreshEval for negate.
        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            return (-1.0)*child_->Eval<dbl>(diff_variable);
        }
        
        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            return -child_->Eval<mpfr>(diff_variable);
        }

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
					if (*(children_mult_or_div_.begin() + (iter-children_.begin())+1)) { // TODO i think this +1 is wrong... dab
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
        virtual std::shared_ptr<Node> Differentiate() override;



         /**
		Compute the degree of a node.  For trig functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return nan.
        */
		virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		{
			int deg = 0;
			for (auto iter = children_.begin(); iter!= children_.end(); iter++)
			{
				// if the child node is a differential coming from another variable, then this degree should be 0, end of story.

				auto factor_deg = (*iter)->Degree(v);

				auto is_it_a_differential = std::dynamic_pointer_cast<Differential>(*iter);
				if (is_it_a_differential)
					if (is_it_a_differential->GetVariable()!=v)
					{
						// std::cout << *this << " deg " << 0 << "\n";
						return 0;
					}
					
				

				

				if (factor_deg<0)
					return factor_deg;
				else if (factor_deg!=0 && !*(children_mult_or_div_.begin() + (iter-children_.begin()) ) )
					return -1;
				else
					deg+=factor_deg;
			}
			// std::cout << *this << " deg " << deg << "\n";
			return deg;
		}

		virtual void Homogenize(std::vector< std::shared_ptr< Variable > > const& vars, std::shared_ptr<Variable> const& homvar) override
		{
			for (auto iter: children_)
			{
				iter->Homogenize(vars, homvar);
			}
		}

	protected:

		// Specific implementation of FreshEval for mult and divide.
		//  If child_mult_ = true, then multiply, else divide
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




		virtual void Homogenize(std::vector< std::shared_ptr< Variable > > const& vars, std::shared_ptr<Variable> const& homvar) override
		{
			if (exponent_->Degree(vars)==0)
			{
				base_->Homogenize(vars, homvar);
			}
			else{
				throw std::runtime_error("asking for homogenization on non-polynomial node");
				//TODO: this will leave the system in a broken state, partially homogenized...
			}
		}


		virtual ~PowerOperator() = default;










	protected:
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


	
	// begin the overload of operators


	inline std::shared_ptr<bertini::Node> pow(const std::shared_ptr<bertini::Node> & N, const std::shared_ptr<bertini::Node> & p)
	{
		return std::make_shared<bertini::PowerOperator>(N,p);
	}

	inline std::shared_ptr<bertini::Node> pow(const std::shared_ptr<bertini::Node> & N, double p)
	{
		return std::make_shared<bertini::PowerOperator>(N,std::make_shared<bertini::Number>(p));
	}


	// Node -> UnaryOperator -> SqrtOperator
	// Description: This class represents the square root function.  FreshEval method
	// is defined for square root and takes the square root of the child node.
	class SqrtOperator : public  virtual UnaryOperator
	{
	public:
		
		SqrtOperator(){}
		
		SqrtOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
		

		
		virtual void print(std::ostream & target) const override
		{
			target << "sqrt(";
			child_->print(target);
			target << ")";
		}
        
        
        /**
         Differentiates the square root function.
         */
        virtual std::shared_ptr<Node> Differentiate() override
        {
            auto ret_mult = std::make_shared<MultOperator>();
            ret_mult->AddChild(std::make_shared<PowerOperator>(child_, std::make_shared<Number>(-0.5)));
            ret_mult->AddChild(child_->Differentiate());
            ret_mult->AddChild(std::make_shared<Number>(0.5));
            return ret_mult;
        }



		/**
		Compute the degree with respect to a single variable.

		For the square root function, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
	    */
	    virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
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




		virtual ~SqrtOperator() = default;
		
	protected:
		// Specific implementation of FreshEval for negate.
//		dbl FreshEval(dbl) override
//		{
//			return sqrt(child_->Eval<dbl>());
//		}
//		
//		mpfr FreshEval(mpfr) override
//		{
//			return sqrt(child_->Eval<mpfr>());
//		}
        
        // Specific implementation of FreshEval for negate.
        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            return sqrt(child_->Eval<dbl>(diff_variable));
        }
        
        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            return sqrt(child_->Eval<mpfr>(diff_variable));
        }

	};



	inline std::shared_ptr<bertini::Node> sqrt(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::SqrtOperator>(N);
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


	// Node -> UnaryOperator -> ExpOperator
	// Description: This class represents the exponential function.  FreshEval method
	// is defined for exponential and takes the exponential of the child node.
	class ExpOperator : public  virtual UnaryOperator
	{
	public:
		
		ExpOperator(){}
		
		ExpOperator(const std::shared_ptr<Node> & N) : UnaryOperator(N)
		{};
	 

		
		
		virtual void print(std::ostream & target) const override
		{
			target << "exp(";
			child_->print(target);
			target << ")";
		}
        
        
        /**
         Differentiates the exponential function.
         */
        virtual std::shared_ptr<Node> Differentiate() override
        {
            auto ret_mult = std::make_shared<MultOperator>();
            ret_mult->AddChild(std::make_shared<ExpOperator>(child_));
            ret_mult->AddChild(child_->Differentiate());
            return ret_mult;
        }



		/**
		Compute the degree with respect to a single variable.

		For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
	    */
	    virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
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

		
		virtual ~ExpOperator() = default;
		
	protected:
        
        // Specific implementation of FreshEval for exponentiate.
        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            return exp(child_->Eval<dbl>(diff_variable));
        }
        
        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            return exp(child_->Eval<mpfr>(diff_variable));
        }

	};


	inline std::shared_ptr<bertini::Node> exp(const std::shared_ptr<bertini::Node> & N)
	{
		return std::make_shared<bertini::ExpOperator>(N);
	}


} // re: namespace bertini






#endif
