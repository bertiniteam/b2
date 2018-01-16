//This file is part of Bertini 2.
//
//python/function_tree.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/function_tree.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/function_tree.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//  Danielle Brake
//  UWEC
//  Spring 2018
//

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

