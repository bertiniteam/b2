//This file is part of Bertini 2.
//
//python/src/minieigen_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/src/minieigen_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/src/minieigen_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  University of Wisconsin - Eau Claire
//  Spring 2018
//

#include "minieigen_export.hpp"

namespace bertini{
	namespace python{

void ExportMinieigen()
{

	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".minieigen");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("minieigen") = new_submodule;

	scope new_submodule_scope = new_submodule;
	new_submodule_scope.attr("__doc__") = "A home for Eigen functionality exposed to PyBertini through Minieigen";


	// minieigen methods for converting python sequence into Eigen vector or matrix
	custom_VectorAnyAny_from_sequence<Eigen::Matrix<dbl,Eigen::Dynamic,1>>();
	custom_VectorAnyAny_from_sequence<Eigen::Matrix<mpfr,Eigen::Dynamic,1>>();
	custom_MatrixAnyAny_from_sequence<Eigen::Matrix<dbl,Eigen::Dynamic,Eigen::Dynamic>>();
	custom_MatrixAnyAny_from_sequence<Eigen::Matrix<mpfr,Eigen::Dynamic,Eigen::Dynamic>>();
	
	
	// Eigen Vector of type dbl or mpfr
	class_<Eigen::Matrix<dbl,Eigen::Dynamic,1>>(
			"VectorXd",
			"Double precision complex vector of arbitrary runtime size, built on MiniEigen, on top of Eigen",
			py::init<>())
		.def(py::init<int>())
		.def(VectorVisitor<Eigen::Matrix<dbl,Eigen::Dynamic,1>>());

	class_<Eigen::Matrix<mpfr,Eigen::Dynamic,1>>(
			"VectorXmp",
			"Runtime-adjustable multiple precision complex vector of arbitrary runtime size, built on MiniEigen, on top of Eigen",
			py::init<>())
		.def(py::init<int>())
		.def(VectorVisitor<Eigen::Matrix<mpfr,Eigen::Dynamic,1>>());
	
	// Eigen Matrix of type dbl or mpfr
	class_<Eigen::Matrix<dbl,Eigen::Dynamic,Eigen::Dynamic>>(
			"MatrixXd",
			"Double precision complex matrix type, runtime-adjustable size",
			py::init<>())
		.def(py::init<int,int>())
		.def(MatrixVisitor<Eigen::Matrix<dbl,Eigen::Dynamic,Eigen::Dynamic>>());

	class_<Eigen::Matrix<mpfr,Eigen::Dynamic,Eigen::Dynamic>>(
			"MatrixXmp",
			"Variable precision complex matrix of runtime-adjustable size",
			py::init<>())
		.def(py::init<int,int>())
		.def(MatrixVisitor<Eigen::Matrix<mpfr,Eigen::Dynamic,Eigen::Dynamic>>());

};

}} // namespaces