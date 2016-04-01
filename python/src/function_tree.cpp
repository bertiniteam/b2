
#include "function_tree.hpp"

using namespace boost::python;

using dbl = std::complex<double>;
using mpfr = bertini::complex;

// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NodeEvalOverloadsMpfr, bertini::node::Node::template Eval<mpfr>, 0, 1);
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(SumAddChildOverloads, bertini::node::SumOperator::AddChild, 1, 2);
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MultAddChildOverloads, bertini::node::MultOperator::AddChild, 1, 2);



namespace bertini{
	namespace python{
		using Nodeptr = std::shared_ptr<node::Node>;
		

		void SetupFunctionTree()
		{
			// Tell Python that pointers to derived Nodes can be used as Node pointers
			implicitly_convertible<std::shared_ptr<node::Float>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::special_number::Pi>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::special_number::E>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::Variable>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::Differential>, Nodeptr>();
			
			implicitly_convertible<std::shared_ptr<node::SumOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::MultOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::PowerOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::NegateOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::IntegerPowerOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::SqrtOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::ExpOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::LogOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::TrigOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::SinOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::CosOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::TanOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::ArcSinOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::ArcCosOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::ArcTanOperator>, Nodeptr>();
			
			implicitly_convertible<std::shared_ptr<node::Function>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<node::Jacobian>, Nodeptr>();
			
			
			
			// Expose the deque containers
			class_< std::deque< std::shared_ptr< node::Variable > > >("VariableGroup")
			.def(vector_indexing_suite< std::deque< std::shared_ptr< node::Variable > >, true >())
			;
		}
		
		
	}
}

