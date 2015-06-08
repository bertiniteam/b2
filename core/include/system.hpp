//This file is part of Bertini 2.0.
//
//system.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//system.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with system.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015
//
// system.hpp:  provides the bertini::system class.




#include <vector>
#include "function_tree.hpp"
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/number.hpp>



namespace bertini {
	
	
	/**
	 The fundamental polynomial system class for Bertini2.
	 */
	class System{
		
		using Point = std::vector<bertini::complex>; // NOT CORRECT.  should be using Eigen3 vectors, probably
		
		
		using Fn = std::shared_ptr<Function>;
		using Var = std::shared_ptr<Variable>;
	public:
		
		
	private:
		
		
		std::vector< Fn > functions_;
		std::vector< Fn > subfunctions_;
		std::vector< Fn > explicit_parameters_;
		
		
		std::vector< Var > variables_;
		std::vector< Var > implicit_paramters_;
		
		
		Var path_variable_;
		
		std::vector< Fn > constants_;
		
		unsigned precision_;
		
		Point solutions_;
	};
	
}