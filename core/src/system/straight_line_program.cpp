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
//
// written with the help of UWEC student Mike Mumm, summer 2021

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
		if (new_precision==DoublePrecision()){
			// nothing
		}
		else{
			auto& mem = std::get<std::vector<mpfr_complex>>(memory_);



			for (auto& p: true_values_of_numbers_)
			{
				auto& n = std::get<Nd>(p);
				auto& loc = std::get<size_t>(p);

				mem[loc] = n->Eval<mpfr_complex>();
			}

			for (auto& n : mem)
				Precision(n, new_precision);
		}


	}




	void StraightLineProgram::AddInstruction(Operation binary_op, size_t in_loc1, size_t in_loc2, size_t out_loc){
		this->instructions_.push_back(binary_op);
		this->instructions_.push_back(in_loc1);
		this->instructions_.push_back(in_loc2);
		this->instructions_.push_back(out_loc);
	}



	void StraightLineProgram::AddInstruction(Operation unary_op, size_t in_loc, size_t out_loc){
		this->instructions_.push_back(unary_op);
		this->instructions_.push_back(in_loc);
		this->instructions_.push_back(out_loc);
	}

	void StraightLineProgram::AddNumber(Nd const num, size_t loc){
		this->true_values_of_numbers_.push_back(std::pair<Nd,size_t>(num, loc));
	}



}







// stuff for SLPCompiler
namespace bertini{
	using SLP = StraightLineProgram;


	void SLPCompiler::Visit(node::Variable const& n){
		std::cout << "Variables were added  to memory at the beginning of compilation" << std::endl;
	}


	// wtb: factor out this pattern
	void SLPCompiler::Visit(node::Integer const& n){
		this->DealWithNumber(n); // that sweet template magic.  see slp.hpp for the definition of this template function
	}

	void SLPCompiler::Visit(node::Float const& n){
		this->DealWithNumber(n);
	}

	void SLPCompiler::Visit(node::Rational const& n){
		this->DealWithNumber(n);
	}











	void SLPCompiler::Visit(node::Function const & f){
		// put the location of the accepted node into memory, and copy into an output location.
		auto& n = f.entry_node();
		size_t location_entry;

		if (this->locations_encountered_symbols_.find(n) == this->locations_encountered_symbols_.end())
			n->Accept(*this); // think of calling Compile(n)

		location_entry = this->locations_encountered_symbols_[n];

		size_t location_this_node = locations_encountered_symbols_[std::make_shared<node::Function>(f)];

		slp_under_construction_.AddInstruction(Assign, location_entry, location_this_node);
	}


	// arithmetic
	void SLPCompiler::Visit(node::SumOperator const & op){
		std::cout << "visiting SumOperator: " << std::endl;

		// this loop
		// gets the locations of all the things we're going to add up.
		std::vector<size_t> operand_locations;
		for (auto& n : op.Operands()){

			if (this->locations_encountered_symbols_.find(n)!=this->locations_encountered_symbols_.end())
				n->Accept(*this); // think of calling Compile(n)

			operand_locations.push_back(this->locations_encountered_symbols_[n]);
		}

		const auto& signs = op.GetSigns();

		// seed the loop by adding together the first two things, which must exist by hypothesis
		assert(op.Operands().size() >=2);

		size_t prev_result_loc = next_available_mem_;

		if (signs[0])
			slp_under_construction_.AddInstruction(Add,operand_locations[0],operand_locations[1], next_available_mem_++);
		else
			slp_under_construction_.AddInstruction(Subtract,operand_locations[0],operand_locations[1], next_available_mem_++);

		// this loop
		// actually does the additions for the rest of the operands

		for (size_t ii{2}; ii<op.Operands().size(); ++ii){
			if (signs[ii-1])
				slp_under_construction_.AddInstruction(Add,operand_locations[ii],prev_result_loc,next_available_mem_++);
			else
				slp_under_construction_.AddInstruction(Subtract,operand_locations[ii],prev_result_loc,next_available_mem_++);
		}
		this->locations_encountered_symbols_[std::make_shared<node::SumOperator>(op)] =  next_available_mem_ - 1; //the loop before this made the current available memory available so that is why  we have the -1

		// for each pair, look up the node in the memory structure in SLP



			// add a SUM op to the instructions.

			// add/subtract the nodes together.

			// naive method: just loop over all operands, add/sub them up.

			// improved option?: we could do this in a way to minimize numerical error.  the obvious loop is not good for accumulation of error.  instead, use Pairwise summation
			// do pairs (0,1), (2,3), etc,
			// *then* add those temp vals together, until get to end.
			// see https://en.wikipedia.org/wiki/Pairwise_summation


		// loop over pairs of operands, not one by one as is currently done.
	}









	void SLPCompiler::Visit(node::Jacobian const& n){
		std::cout << "unimplemented visit to node of type Jacobian" << std::endl;


	}

	void SLPCompiler::Visit(node::Differential const& n){
		std::cout << "unimplemented visit to node of type Differential" << std::endl;
	}



	void SLPCompiler::Visit(node::MultOperator const & op){
		std::cout << "visiting MultOperator: " << std::endl;
		std::vector<size_t> operand_locations;
		for (auto& n : op.Operands()){

			if (this->locations_encountered_symbols_.find(n)!=this->locations_encountered_symbols_.end())
				n->Accept(*this); // think of calling Compile(n)

			operand_locations.push_back(this->locations_encountered_symbols_[n]);
		}
		// multiply/divide the nodes together.  naive method is fine.
		const auto& MultOrDiv = op.GetMultOrDiv();// true is multiply and false is divide

		assert(op.Operands().size() >=2);

		size_t prev_result_loc = next_available_mem_;

		if (MultOrDiv[0])
			slp_under_construction_.AddInstruction(Multiply,operand_locations[0],operand_locations[1], next_available_mem_++);
		else
			slp_under_construction_.AddInstruction(Divide,operand_locations[0],operand_locations[1], next_available_mem_++);

		// this loop
		// actually does the additions for the rest of the operands

		for (size_t ii{2}; ii<op.Operands().size(); ++ii){
			if (MultOrDiv[ii-1])
				slp_under_construction_.AddInstruction(Multiply,operand_locations[ii],prev_result_loc,next_available_mem_++);
			else
				slp_under_construction_.AddInstruction(Divide,operand_locations[ii],prev_result_loc,next_available_mem_++);
		}
		this->locations_encountered_symbols_[std::make_shared<node::MultOperator>(op)] =  next_available_mem_ - 1; //the loop before this made the current available memory available so that is why  we have the -1

	}

	void SLPCompiler::Visit(node::IntegerPowerOperator const& n){
		std::cout << "unimplemented visit to node of type IntegerPowerOperator" << std::endl;
		auto expo = n.exponent(); //integer

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		auto e = std::make_shared<node::Integer>(expo);
		this->DealWithNumber(*e);
		auto location_exponent  = locations_encountered_symbols_[e];


		this->locations_encountered_symbols_[std::make_shared<node::IntegerPowerOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Power,location_operand,location_exponent, next_available_mem_++);
		//no node, just integer.
	}

	void SLPCompiler::Visit(node::PowerOperator const& n){
		//get location of base and power then add instruction

		const auto& base = n.GetBase();
		const auto& exponent = n.GetExponent();

		if (this->locations_encountered_symbols_.find(base) == this->locations_encountered_symbols_.end())
			base->Accept(*this); // think of calling Compile(n)

		if (this->locations_encountered_symbols_.find(exponent) == this->locations_encountered_symbols_.end())
			exponent->Accept(*this); // think of calling Compile(n)

		auto loc_base = locations_encountered_symbols_[base];
		auto loc_exponent = locations_encountered_symbols_[exponent];

		this->locations_encountered_symbols_[std::make_shared<node::PowerOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Power, loc_base, loc_exponent, next_available_mem_++);



	}

	void SLPCompiler::Visit(node::ExpOperator const& n){

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[std::make_shared<node::ExpOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Exp,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::LogOperator const& n){

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[std::make_shared<node::LogOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Log,location_operand, next_available_mem_++);
	}


	// the trig operators
	void SLPCompiler::Visit(node::SinOperator const& n){

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[std::make_shared<node::SinOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Sin,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::ArcSinOperator const& n){

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[std::make_shared<node::ArcSinOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Asin,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::CosOperator const& n){

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[std::make_shared<node::CosOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Cos,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::ArcCosOperator const& n){

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[std::make_shared<node::ArcCosOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Acos,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::TanOperator const& n){

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[std::make_shared<node::TanOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Tan,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::ArcTanOperator const& n){

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[std::make_shared<node::ArcTanOperator>(n)] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Atan,location_operand, next_available_mem_++);
	}



	SLP SLPCompiler::Compile(System const& sys){
		this->Clear();

		std::cout << "Compiling system" << std::endl;


			// deal with variables

		std::cout << "dealing with variables" << std::endl;

			// 1. ADD VARIABLES
		auto variable_groups = sys.VariableGroups();
		auto variable_group_sizes = sys.VariableGroupSizes();

		unsigned variable_counter{0};
		for (int ii=0; ii<variable_group_sizes.size(); ++ii){
			auto vg = variable_groups[ii];
			auto s = variable_group_sizes[ii];

			for (int jj=0; jj<s; ++jj){
				locations_encountered_symbols_[ vg[jj] ] = next_available_mem_++;
				variable_counter += s;
			}
		}

		slp_under_construction_.num_variables_ = variable_counter;

			// deal with path variable
		if (sys.HavePathVariable())
		{
				// do this action only if the system has a path variable defined
			locations_encountered_symbols_[ sys.GetPathVariable() ] = next_available_mem_++;
			slp_under_construction_.has_path_variable_ = true;
		}


		std::cout << "making space in memory for functions" << std::endl;
			// make space for functions and derivatives
			// 3. ADD FUNCTIONS
		for  (int ii = 0; ii < sys.NumTotalFunctions(); ++ii) {
			auto f = sys.Function(ii);
			locations_encountered_symbols_[f] = next_available_mem_++;
		}


		
		std::cout << "making space in memory for derivatives" << std::endl;
			// always do derivatives with respect to space variables
			// 4. ADD SPACE VARIABLE DERIVATIVES

		if (sys.HavePathVariable())
		{
				// we need derivatives with respect to time only if the system has a path variable defined
				// 5. ADD TIME VARIABLE DERIVATIVES
		}

		std::cout << "visiting functions" << std::endl;
		for (int ii = 0; ii < sys.NumTotalFunctions(); ++ii)
		{
			auto f = sys.Function(ii);

			std::cout << *(f) << std::endl;
			f->Accept(*this);


			std::cout << "post visit function" << std::endl;
			/* code */
		}

		return slp_under_construction_;
	}

	void SLPCompiler::Clear(){
		next_available_mem_ = 0;
		locations_encountered_symbols_.clear();
		slp_under_construction_ = SLP();
	}
}
