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
// michael mumm, university of wisconsin eau claire

#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/system.hpp"



// SLP Stuff
namespace bertini{


	// the constructor
	StraightLineProgram::StraightLineProgram(System const& sys){
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
		std::cout << "added binary instruction, op " << binary_op << " args: " << in_loc1 << ", " << in_loc2 << " dest: " << out_loc << std::endl;
	}



	void StraightLineProgram::AddInstruction(Operation unary_op, size_t in_loc, size_t out_loc){
		this->instructions_.push_back(unary_op);
		this->instructions_.push_back(in_loc);
		this->instructions_.push_back(out_loc);
		std::cout << "added unary instruction, op " << unary_op << " arg: " << in_loc << " dest: " << out_loc <<  std::endl;

	}

	void StraightLineProgram::AddNumber(Nd const num, size_t loc){
		this->true_values_of_numbers_.push_back(std::pair<Nd,size_t>(num, loc));
		std::cout << "added number " << *num << " at " << loc << std::endl;
	}



	std::ostream& operator <<(std::ostream& out, const StraightLineProgram & s){
		out << "\n\n#fns: " << s.NumFunctions() << " #vars: " << s.NumVariables() << std::endl;

		out << "instructions: " << std::endl;
		for (auto i: s.instructions_)
			out << i << " ";
		out << std::endl;

		out << "numbers: " << std::endl;
		for (auto const& x : s.true_values_of_numbers_)
		    out << *(x.first)  << ':' << x.second << std::endl;
		out << std::endl << std::endl;

		return out;
	}


	template<typename NumT>
	void StraightLineProgram::CopyNumbersIntoMemory() const
	{
		for (auto const& x: true_values_of_numbers_){
			GetMemory<NumT>()[x.second] = (x.first)->Eval<NumT>();
			if (std::is_same<NumT,mpfr_complex>::value)
				Precision(GetMemory<NumT>()[x.second], this->precision_);
		}
	}

	template void StraightLineProgram::CopyNumbersIntoMemory<dbl_complex>() const;
	template void StraightLineProgram::CopyNumbersIntoMemory<mpfr_complex>() const;

}







// stuff for SLPCompiler
namespace bertini{
	using SLP = StraightLineProgram;


	void SLPCompiler::Visit(node::Variable const& n){
		std::cout << "visiting variable " << n << " at location " << locations_encountered_symbols_[n.shared_from_this()] << std::endl;
	}


	// wtb: factor out this pattern
	void SLPCompiler::Visit(node::Integer const& n){
		std::cout << "visiting Integer " << n << std::endl;
		this->DealWithNumber(n); // that sweet template magic.  see slp.hpp for the definition of this template function
	}

	void SLPCompiler::Visit(node::Float const& n){
		std::cout << "visiting Float " << n << std::endl;
		this->DealWithNumber(n);
	}

	void SLPCompiler::Visit(node::Rational const& n){
		std::cout << "visiting Integer " << n << std::endl;
		this->DealWithNumber(n);
	}





	void SLPCompiler::Visit(node::Jacobian const& n){
		throw std::runtime_error("unimplemented visit to node of type Jacobian");


	}

	void SLPCompiler::Visit(node::Differential const& n){
		throw std::runtime_error("unimplemented visit to node of type Differential");
	}





	void SLPCompiler::Visit(node::Function const & f){
		std::cout << "visiting function " << f << std::endl;
		// put the location of the accepted node into memory, and copy into an output location.
		auto& n = f.entry_node();
		size_t location_entry;

		if (this->locations_encountered_symbols_.find(n) == this->locations_encountered_symbols_.end())
			n->Accept(*this); // think of calling Compile(n)

		location_entry = this->locations_encountered_symbols_[n];

		size_t location_this_node;
		if (this->locations_encountered_symbols_.find(f.shared_from_this()) == this->locations_encountered_symbols_.end()){
			location_this_node = next_available_mem_;
			locations_encountered_symbols_[f.shared_from_this()] = next_available_mem_++;
		}
		else
			location_this_node = locations_encountered_symbols_[f.shared_from_this()];


		slp_under_construction_.AddInstruction(Assign, location_entry, location_this_node);
		std::cout << "added function node, " << f << ", copying from " << location_entry << " to " << location_this_node << std::endl;
	}


	// arithmetic
	void SLPCompiler::Visit(node::SumOperator const & n){
		std::cout << "visiting SumOperator: " << n << std::endl;

		// this loop
		// gets the locations of all the things we're going to add up.
		std::vector<size_t> operand_locations;
		for (auto& n : n.Operands()){

			if (this->locations_encountered_symbols_.find(n)==this->locations_encountered_symbols_.end())
				n->Accept(*this); // think of calling Compile(n)

			operand_locations.push_back(this->locations_encountered_symbols_[n]);
			std::cout << "sum operand is: " << n << std::endl;
		}

		  
		const auto& signs = n.GetSigns();
		size_t prev_result_loc; // for tracking where the output of the previous iteration went

		// seed the loop.
		if (signs[0])
			prev_result_loc = operand_locations[0];
		else{
			slp_under_construction_.AddInstruction(Negate, operand_locations[0], next_available_mem_);
			prev_result_loc = next_available_mem_++;
		}
		

		// this loop
		// does the additions for the rest of the operands
		for (size_t ii{1}; ii<n.Operands().size(); ++ii){
			if (signs[ii])
				slp_under_construction_.AddInstruction(Add,prev_result_loc,operand_locations[ii],next_available_mem_);
			else
				slp_under_construction_.AddInstruction(Subtract,prev_result_loc,operand_locations[ii],next_available_mem_);

			prev_result_loc = next_available_mem_++;
		}

		this->locations_encountered_symbols_[n.shared_from_this()] =  prev_result_loc;

			// improved option?: we could do this in a way to minimize numerical error.  the obvious loop is not good for accumulation of error.  instead, use Pairwise summation
			// do pairs (0,1), (2,3), etc,
			// *then* add those temp vals together, until get to end.
			// see https://en.wikipedia.org/wiki/Pairwise_summation
	}







	void SLPCompiler::Visit(node::MultOperator const & n){
		std::cout << "visiting MultOperator: " << n << std::endl;


		// this loop
		// gets the locations of all the things we're going to add up.
		std::vector<size_t> operand_locations;
		for (auto& n : n.Operands()){

			if (this->locations_encountered_symbols_.find(n)==this->locations_encountered_symbols_.end())
				n->Accept(*this); // think of calling Compile(n)

			operand_locations.push_back(this->locations_encountered_symbols_[n]);
			std::cout << "mult operand is: " << n << std::endl;
		}

		  
		const auto& mult_or_div = n.GetMultOrDiv();// true is multiply and false is divide

		size_t prev_result_loc; // for tracking where the output of the previous iteration went

		// seed the loop.
		if (mult_or_div[0])
			prev_result_loc = operand_locations[0];
		else{ 
			// this case is reciprocation of the first operand

			// this code sucks.  really, there should be a bank of integers that we pull from, instead of many copies of the same integer.
			auto one = std::make_shared<node::Integer>(1);
			this->DealWithNumber(*one);
			auto location_one  = locations_encountered_symbols_[one];

			slp_under_construction_.AddInstruction(Divide, location_one, operand_locations[0], next_available_mem_);
			prev_result_loc = next_available_mem_++;
		}
		

		// this loop
		// does the additions for the rest of the operands
		for (size_t ii{1}; ii<n.Operands().size(); ++ii){
			if (mult_or_div[ii])
				slp_under_construction_.AddInstruction(Multiply,prev_result_loc,operand_locations[ii],next_available_mem_);
			else
				slp_under_construction_.AddInstruction(Divide,prev_result_loc,operand_locations[ii],next_available_mem_);

			prev_result_loc = next_available_mem_++;
		}

		this->locations_encountered_symbols_[n.shared_from_this()] =  prev_result_loc;

	}





	void SLPCompiler::Visit(node::IntegerPowerOperator const& n){
		std::cout << "visiting IntegerPowerOperator" << n << std::endl;
		auto expo = n.exponent(); //integer

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)


		// this code sucks, because it results in many copies of the same few integers over and over.  
		// instead, we should add a bank of integers, and only add a new one if needed.
		auto location_operand = locations_encountered_symbols_[operand];
		auto e = std::make_shared<node::Integer>(expo);
		this->DealWithNumber(*e);
		auto location_exponent  = locations_encountered_symbols_[e];


		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Power,location_operand,location_exponent, next_available_mem_++);
		//no node, just integer.
	}


	void SLPCompiler::Visit(node::PowerOperator const& n){
		std::cout << "visiting PowerOperator" << n << std::endl;
		//get location of base and power then add instruction

		const auto& base = n.GetBase();
		const auto& exponent = n.GetExponent();

		if (this->locations_encountered_symbols_.find(base) == this->locations_encountered_symbols_.end())
			base->Accept(*this); // think of calling Compile(n)

		if (this->locations_encountered_symbols_.find(exponent) == this->locations_encountered_symbols_.end())
			exponent->Accept(*this); // think of calling Compile(n)

		auto loc_base = locations_encountered_symbols_[base];
		auto loc_exponent = locations_encountered_symbols_[exponent];

		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Power, loc_base, loc_exponent, next_available_mem_++);



	}

	void SLPCompiler::Visit(node::ExpOperator const& n){
		std::cout << "visiting ExpOperator" << n << std::endl;

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Exp,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::LogOperator const& n){
		std::cout << "visiting LogOperator" << n << std::endl;

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Log,location_operand, next_available_mem_++);
	}


	// the trig operators
	void SLPCompiler::Visit(node::SinOperator const& n){
		std::cout << "visiting SinOperator" << n << std::endl;

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Sin,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::ArcSinOperator const& n){
		std::cout << "visiting ArcSinOperator" << n << std::endl;

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Asin,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::CosOperator const& n){
		std::cout << "visiting CosOperator" << n << std::endl;

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Cos,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::ArcCosOperator const& n){
		std::cout << "visiting ArcCosOperator" << n << std::endl;

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Acos,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::TanOperator const& n){
		std::cout << "visiting TanOperator" << n << std::endl;

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Tan,location_operand, next_available_mem_++);
	}

	void SLPCompiler::Visit(node::ArcTanOperator const& n){
		std::cout << "visiting ArcTanOperator" << n << std::endl;

		auto operand = n.Operand();
		if (this->locations_encountered_symbols_.find(operand) == this->locations_encountered_symbols_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_symbols_[operand];
		this->locations_encountered_symbols_[n.shared_from_this()] =  next_available_mem_;
		slp_under_construction_.AddInstruction(Atan,location_operand, next_available_mem_++);
	}



	SLP SLPCompiler::Compile(System const& sys){
		this->Clear();

		std::cout << "Compiling system" << sys << std::endl;


		this->slp_under_construction_.precision_ = DefaultPrecision();

		// deal with variables
		std::cout << "adding variables to memory" << std::endl;

			// 1. ADD VARIABLES
		auto variable_groups = sys.VariableGroups();
		auto variable_group_sizes = sys.VariableGroupSizes();
		slp_under_construction_.input_locations_.Variables = next_available_mem_;
		unsigned variable_counter{0};
		for (int ii=0; ii<variable_group_sizes.size(); ++ii){
			auto vg = variable_groups[ii];
			auto s = variable_group_sizes[ii];

			for (int jj=0; jj<s; ++jj){
				locations_encountered_symbols_[ vg[jj] ] = next_available_mem_++;
				variable_counter += 1;
			}
		}
		slp_under_construction_.number_of_.Variables = variable_counter;
		slp_under_construction_.input_locations_.Variables = variable_counter;

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
		slp_under_construction_.number_of_.Functions = sys.NumFunctions();
		slp_under_construction_.output_locations_.Functions = next_available_mem_;
		for (auto f: sys.GetFunctions())
			locations_encountered_symbols_[f] = next_available_mem_++;





		std::cout << "making space in memory for space derivatives" << std::endl;
		// always have space derivatives

		auto ds_dx = sys.GetSpaceDerivatives(); // a linear object, so can just run down the object
		slp_under_construction_.number_of_.Jacobian = ds_dx.size();
		slp_under_construction_.output_locations_.Jacobian = next_available_mem_;
		for (auto n: ds_dx)
			locations_encountered_symbols_[n] = next_available_mem_++;


		// sometimes have time derivatives
		if (sys.HavePathVariable()) {
			std::cout << "making space in memory for time derivatives" << std::endl;
			auto ds_dt = sys.GetTimeDerivatives();  // a linear object, so can just run down the object
			slp_under_construction_.number_of_.TimeDeriv = ds_dt.size();
			slp_under_construction_.output_locations_.TimeDeriv = next_available_mem_;
			for (auto n: ds_dt)
				locations_encountered_symbols_[n] = next_available_mem_++;
		}






		std::cout << "visiting functions" << std::endl;
		for (auto f: sys.GetFunctions())
		{
			f->Accept(*this);

			// post visit function
			/* code */
		}



		std::cout << "visiting space derivatives" << std::endl;
		// always do derivatives with respect to space variables
		// 4. ADD SPACE VARIABLE DERIVATIVES
		for (auto n: ds_dx)
			n->Accept(*this);




		// sometimes have time derivatives
		if (sys.HavePathVariable()) {
			std::cout << "visiting time derivatives" << std::endl;
			// we need derivatives with respect to time only if the system has a path variable defined
			// 5. ADD TIME VARIABLE DERIVATIVES

			auto ds_dt = sys.GetTimeDerivatives();  // a linear object, so can just run down the object
			for (auto n: ds_dt)
				n->Accept(*this);
		}

		slp_under_construction_.GetMemory<dbl_complex>().resize(next_available_mem_);
		slp_under_construction_.GetMemory<mpfr_complex>().resize(next_available_mem_);

		slp_under_construction_.CopyNumbersIntoMemory<dbl_complex>();
		slp_under_construction_.CopyNumbersIntoMemory<mpfr_complex>();

		return slp_under_construction_;
	}

	void SLPCompiler::Clear(){
		next_available_mem_ = 0;
		locations_encountered_symbols_.clear();
		slp_under_construction_ = SLP();
	}
}
