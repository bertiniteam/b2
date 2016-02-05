//
//  root_export.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/5/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_root_export_hpp
#define Xcode_b2_root_export_hpp
#include <bertini2/function_tree/roots/function.hpp>
#include <bertini2/function_tree/roots/jacobian.hpp>

#include "python_common.hpp"

namespace bertini{
	namespace python{
		
		using namespace boost::python;
		using namespace bertini::node;
		
		using Node = Node;
		using Nodeptr = std::shared_ptr<Node>;
		
		void ExportRoots();
		
		template<typename NodeBaseT>
		class FunctionVisitor: public def_visitor<FunctionVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		};

		
		template<typename NodeBaseT>
		class JacobianVisitor: public def_visitor<JacobianVisitor<NodeBaseT> >
		{
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
		};

		
	}
}

#endif
