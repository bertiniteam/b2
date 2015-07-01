
#include "function_tree/operators/sum_operator.hpp"


namespace bertini{

	using ::pow;

	void SumOperator::Homogenize(std::vector< std::shared_ptr< Variable > > const& vars, std::shared_ptr<Variable> const& homvar)
	{
		
	
	
		// first homogenize each summand.
		for (auto iter: children_)
		{
			iter->Homogenize(vars, homvar);
		}

		// then, homogenize this sum.

		// compute the highest degree among all summands.
		int maxdegree = 0;
		std::vector<int> term_degrees;
		// first homogenize each summand.
		for (auto iter: children_)
		{
			auto local_degree = iter->Degree(vars);
			if (local_degree<0)
				throw std::runtime_error("asking for homogenization on non-polynomial node");
				// TODO: this throw would leave the tree in a partially homogenized state.  this is scary.

			term_degrees.push_back(local_degree);
			maxdegree = std::max(maxdegree, local_degree);
		}

		for (auto iter = children_.begin(); iter!=children_.end(); iter++)
		{
			auto degree_deficiency = maxdegree - *(term_degrees.begin() + (iter-children_.begin()));
			if ( degree_deficiency > 0)
			{
				// hold the child temporarily.
				std::shared_ptr<Node> P = std::make_shared<IntegerPowerOperator>(std::dynamic_pointer_cast<Node>(homvar),degree_deficiency);
				std::shared_ptr<Node> M = std::make_shared<MultOperator>(P,std::dynamic_pointer_cast<Node>(*iter));
				
				
				swap(*iter,M);
			}
		}

	}

}

