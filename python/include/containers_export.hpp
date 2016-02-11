//
//  export_containers.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/10/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_export_containers_hpp
#define Xcode_b2_export_containers_hpp


#include <bertini2/function_tree.hpp>

#include "python_common.hpp"




void ExportContainers()
{
	// Expose the deque containers
	class_< bertini::VariableGroup >("VariableGroup")
	.def(vector_indexing_suite< bertini::VariableGroup, true >())
	;
	
	
	class_< std::vector<int> >("ListInt")
	.def(vector_indexing_suite< std::vector<int> , true >())
	;

	class_< std::vector<bertini::VariableGroup> >("ListVariableGroup")
	.def(vector_indexing_suite< std::vector<bertini::VariableGroup> , true >())
	;

	class_< std::vector<std::shared_ptr< bertini::node::Function > > >("ListFunction")
	.def(vector_indexing_suite< std::vector<std::shared_ptr< bertini::node::Function> > , true >())
	;

	class_< std::vector<std::shared_ptr< bertini::node::Jacobian > > >("ListJacobian")
	.def(vector_indexing_suite< std::vector<std::shared_ptr< bertini::node::Jacobian> > , true >())
	;

};








#endif
