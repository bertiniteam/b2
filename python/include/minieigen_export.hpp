//
//  minieigen_export.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/11/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_minieigen_export_hpp
#define Xcode_b2_minieigen_export_hpp

#include <bertini2/function_tree.hpp>

#include "../minieigen/src/common.hpp"
#include "../minieigen/src/expose.hpp"

#include "python_common.hpp"

namespace bertini{
	namespace python{
		
		
		void ExportMinieigen()
		{
			expose_converters(); // in expose-converters.cpp
			
			
			class_<Eigen::Matrix<dbl,Eigen::Dynamic,1>>("VectorXd","/*TODO*/",
														 py::init<>()).def(VectorVisitor<Eigen::Matrix<dbl,Eigen::Dynamic,1>>());
			class_<Eigen::Matrix<mpfr,Eigen::Dynamic,1>>("VectorXmp","/*TODO*/",
														py::init<>()).def(VectorVisitor<Eigen::Matrix<mpfr,Eigen::Dynamic,1>>());
			
			class_<Eigen::Matrix<dbl,Eigen::Dynamic,Eigen::Dynamic>>("MatrixXd","/*TODO*/",
														py::init<>()).def(VectorVisitor<Eigen::Matrix<dbl,Eigen::Dynamic,Eigen::Dynamic>>());
			class_<Eigen::Matrix<mpfr,Eigen::Dynamic,Eigen::Dynamic>>("MatrixXmp","/*TODO*/",
														 py::init<>()).def(VectorVisitor<Eigen::Matrix<mpfr,Eigen::Dynamic,Eigen::Dynamic>>());

		};
		
		
	}
}


#endif
