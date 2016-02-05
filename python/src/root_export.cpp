//
//  root_export.cpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/5/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#include <stdio.h>
#include "root_export.hpp"


namespace bertini{
	namespace python{
		
		
		
		template<typename NodeBaseT>
		template<class PyClass>
		void FunctionVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.add_property("root", &Function::entry_node, &Function::SetRoot)
			.def("ensure_not_empy", &Function::EnsureNotEmpty)
			;
		}

		
		template<typename NodeBaseT>
		template<class PyClass>
		void JacobianVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("evalJd", &Jacobian::template EvalJ<dbl>)
			.def("evalJmp", &Jacobian::template EvalJ<mpfr>)
			;
		}

		
		
		void ExportRoots()
		{
			// Function class
			class_<Function, bases<NamedSymbol>, std::shared_ptr<Function> >("Function", init<>())
			.def(init<std::string>() )
			.def(init<const Nodeptr&>() )
			
			.def(FunctionVisitor<Function>())
			;

			
			// Jacobian class
			class_<Jacobian, bases<Function>, std::shared_ptr<Jacobian> >("Jacobian", init<>())
			.def(init<const Nodeptr&>() )
			
			.def(JacobianVisitor<Jacobian>())
			;
			
		}

	}
}

