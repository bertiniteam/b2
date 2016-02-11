//
//  minieigen_export.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/11/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_minieigen_export_hpp
#define Xcode_b2_minieigen_export_hpp

#include "../minieigen/src/common.hpp"
#include "../minieigen/src/expose.hpp"

#include "python_common.hpp"

namespace bertini{
	namespace python{
		
		
		void ExportMinieigen()
		{
			expose_converters(); // in expose-converters.cpp
			
			expose_vectors();
			expose_matrices(); // must come after vectors
			expose_complex();

		};
		
		
	}
}


#endif
