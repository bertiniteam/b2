//This file is part of Bertini 2.
//
//python/function_tree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/function_tree.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/function_tree.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//
//  python/function_tree.hpp:  Header file for all node types.

#pragma once
#ifndef BERTINI_PYTHON_FUNCTION_TREE_HPP
#define BERTINI_PYTHON_FUNCTION_TREE_HPP

#include "python_common.hpp"
#include <bertini2/function_tree.hpp>

#include "node_export.hpp"
#include "symbol_export.hpp"
#include "operator_export.hpp"
#include "root_export.hpp"



namespace bertini{
	namespace python{

		using namespace bertini::node;
		using Node = Node;
		using Nodeptr = std::shared_ptr<Node>;


		

		void SetupFunctionTree()
		{
			// Tell Python that pointers to derived Nodes can be used as Node pointers when passed to methods
			implicitly_convertible<std::shared_ptr<Float>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<special_number::Pi>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<special_number::E>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<Integer>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<Rational>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<Variable>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<Differential>, Nodeptr>();
			
			implicitly_convertible<std::shared_ptr<SumOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<MultOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<PowerOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<NegateOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<IntegerPowerOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<SqrtOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<ExpOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<LogOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<TrigOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<SinOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<CosOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<TanOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<ArcSinOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<ArcCosOperator>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<ArcTanOperator>, Nodeptr>();
			
			implicitly_convertible<std::shared_ptr<Function>, Nodeptr>();
			implicitly_convertible<std::shared_ptr<Jacobian>, Nodeptr>();
			
		}
		
		
		
		
		
		
	}
}


#endif
