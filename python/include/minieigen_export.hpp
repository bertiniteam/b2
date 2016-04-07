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
//  python/minieigen_export.hpp:  Header file for exposing the eigen data types needed for bertini.




#ifndef Xcode_b2_minieigen_export_hpp
#define Xcode_b2_minieigen_export_hpp

#include <bertini2/function_tree.hpp>
#include <bertini2/eigen_extensions.hpp>

#include "python_common.hpp"

#include "minieigen/src/common.hpp"
#include "minieigen/src/converters.hpp"
#include "minieigen/src/visitors.hpp"





namespace bertini{
	namespace python{
		
		
		void ExportMinieigen()
		{
			// minieigen methods for converting python sequence into Eigen vector or matrix
			custom_VectorAnyAny_from_sequence<Eigen::Matrix<dbl,Eigen::Dynamic,1>>();
			custom_VectorAnyAny_from_sequence<Eigen::Matrix<mpfr,Eigen::Dynamic,1>>();
			custom_MatrixAnyAny_from_sequence<Eigen::Matrix<dbl,Eigen::Dynamic,Eigen::Dynamic>>();
			custom_MatrixAnyAny_from_sequence<Eigen::Matrix<mpfr,Eigen::Dynamic,Eigen::Dynamic>>();
			
			
			// Eigen Vector of type dbl or mpfr
			class_<Eigen::Matrix<dbl,Eigen::Dynamic,1>>("VectorXd","/*TODO*/",
														 py::init<>()).def(VectorVisitor<Eigen::Matrix<dbl,Eigen::Dynamic,1>>());
			class_<Eigen::Matrix<mpfr,Eigen::Dynamic,1>>("VectorXmp","/*TODO*/",
														py::init<>()).def(VectorVisitor<Eigen::Matrix<mpfr,Eigen::Dynamic,1>>());
			
			// Eigen Matrix of type dbl or mpfr
			class_<Eigen::Matrix<dbl,Eigen::Dynamic,Eigen::Dynamic>>("MatrixXd","/*TODO*/",
														py::init<>()).def(MatrixVisitor<Eigen::Matrix<dbl,Eigen::Dynamic,Eigen::Dynamic>>());
			class_<Eigen::Matrix<mpfr,Eigen::Dynamic,Eigen::Dynamic>>("MatrixXmp","/*TODO*/",
														 py::init<>()).def(MatrixVisitor<Eigen::Matrix<mpfr,Eigen::Dynamic,Eigen::Dynamic>>());

		};
		
		
	}
}


#endif
