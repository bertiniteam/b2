// python/function_tree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// python/function_tree.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/function_tree.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Jeb Collins
//  West Texas A&M University
//  Mathematics
//  Fall 2015
//
//
//  python/system.hpp:  the header file for the python interface for System classes.

#ifndef BERTINI_PYTHON_SYSTEM_HPP
#define BERTINI_PYTHON_SYSTEM_HPP

#include "function_tree.hpp"

#include <Eigen/Dense>
#include <vector>



namespace bertini{
	namespace python{
		
		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
		template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
		using VariableGroup = std::deque< std::shared_ptr<node::Variable> >;
		
		using dbl = std::complex<double>;
		using mpfr = bertini::complex;

		using namespace boost::python;
		

		void (bertini::System::*sysAddFunc1)(std::shared_ptr<node::Function> const&) = &bertini::System::AddFunction;
		void (bertini::System::*sysAddFunc2)(std::shared_ptr<node::Node> const&) = &bertini::System::AddFunction;
		
		template<typename T> Vec<T> (bertini::System::*sysEval1)(const Vec<T> &) = &bertini::System::Eval<T>;
		template<typename T> Vec<T> (bertini::System::*sysEval2)(const Vec<T> &, const T &) = &bertini::System::Eval<T>;
		
		template<typename T> Mat<T> (bertini::System::*sysJac1)(const Vec<T> &) = &bertini::System::Jacobian<T>;
		template<typename T> Mat<T> (bertini::System::*sysJac2)(const Vec<T> &, const T &) = &bertini::System::Jacobian<T>;
				
		std::vector<int> (bertini::System::*sysDeg1)() const = &bertini::System::Degrees;
		std::vector<int> (bertini::System::*sysDeg2)(VariableGroup const&) const = &bertini::System::Degrees;
	} // re: python
} // re: bertini

#endif
