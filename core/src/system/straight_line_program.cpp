//This file is part of Bertini 2.
//
//straight_line_program.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//straight_line_program.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with straight_line_program.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/system.hpp"



// SLP Stuff
namespace bertini{
	

	// the constructor
	StraightLineProgram::StraightLineProgram(System const& sys){
		this->num_total_functions_ = sys.NumTotalFunctions();

		std::cout << sys.NumTotalFunctions() << std::endl;

		SLPCompiler compiler;
		*this = compiler.Compile(sys);
	}

	void StraightLineProgram::precision(unsigned new_precision) const{

	}
}



// stuff for SLPCompiler
namespace bertini{
	using SLP = StraightLineProgram;

	void SLPCompiler::Visit(node::Function const & f){
		std::cout << "visiting Function: " << f << std::endl;
		f.entry_node()->Accept(*this);
	}


    void SLPCompiler::Visit(node::SumOperator const & op){
    	std::cout << "visiting SumOperator: " << std::endl;
    }

    void SLPCompiler::Visit(node::Node const & n){
    	std::cout << "visiting generic Node: " << std::endl;
    }

    SLP SLPCompiler::Compile(System const& sys){
    	std::cout << "Compiling system" << std::endl;

    	std::cout << "visiting functions" << std::endl;

    	for (int ii = 0; ii < sys.NumTotalFunctions(); ++ii)
    	{
    		auto f = sys.Function(ii);

    		std::cout << *(f) << std::endl;
    		f->Accept(*this);

    		locations_encountered_symbols_[f] = 0; // this is obviously wrong

    		std::cout << "post visit function" << std::endl;
    		/* code */
    	}

    	return SLP();
    }
}