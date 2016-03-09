//This file is part of Bertini 2.0.
//
// python/function_tree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//This file is part of Bertini 2.0.
//
// python/bertini_python.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/bertini_python.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/bertini_python.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//
//  python/root_export.cpp:  Source file for exposing root nodes to python.




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
			class_<Function, bases<NamedSymbol>, std::shared_ptr<Function> >("Function", init<std::string>())
			.def(init<const Nodeptr&>() )
			
			.def(FunctionVisitor<Function>())
			;

			
			// Jacobian class
			class_<Jacobian, bases<Function>, std::shared_ptr<Jacobian> >("Jacobian", init<const Nodeptr&>())
			
			.def(JacobianVisitor<Jacobian>())
			;
			
		}

	}
}

