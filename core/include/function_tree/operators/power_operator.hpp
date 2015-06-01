
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
			target << "(" << *base_ << "^" << *exponent_ << ")";
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

