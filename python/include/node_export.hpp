//
//  node_export.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 1/30/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_node_export_hpp
#define Xcode_b2_node_export_hpp

#include "node_visitors.hpp"

//#include <bertini2/function_tree/node.hpp>


namespace bertini{
	namespace python{
		
		

		void ExportNode()
		{
			class_<NodeWrap, boost::noncopyable, Nodeptr >("Node", no_init)
			.def(NodeVisitor<NodeWrap>());
			
		};
		
		
	} //namespace python
} // namespace bertini

#endif
