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

	std::string OpcodeToString(Operation op)
	{
		switch (op){
		case Add: return "Add";
		case Subtract: return "Subtract";
		case Multiply: return "Multiply";
		case Divide: return "Divide";
		case Power: return "Power";
		case Exp: return "Exp";
		case Log: return "Log";
		case Negate: return "Negate";
		case Sqrt: return "Sqrt";
		case Sin: return "Sin";
		case Cos: return "Cos";
		case Tan: return "Tan";
		case Asin: return "Asin";
		case Acos: return "Acos";
		case Atan: return "Atan";
		case Assign: return "Assign";
		case IntPower: return "IntPower";
		}
	}


	// the constructor
	StraightLineProgram::StraightLineProgram(System const& sys){
		SLPCompiler compiler;

		*this = compiler.Compile(sys);
	}

	void StraightLineProgram::precision(unsigned new_precision) const{

		if (new_precision==this->precision_){
			return;
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

			this->precision_ = new_precision;
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



	std::ostream& operator <<(std::ostream& out, const StraightLineProgram & s){
		out << "\n\n#fns: " << s.NumFunctions() << " #vars: " << s.NumVariables() << std::endl;
		out << "have path variable: " << s.HavePathVariable() << std::endl;

		out << std::endl << "numbers of things:" << std::endl;
		out << "Functions: " << s.number_of_.Functions << std::endl;
		out << "Variables: " << s.number_of_.Variables << std::endl;
		out << "Jacobian: " << s.number_of_.Jacobian << std::endl;

		if (s.HavePathVariable())
			out << "TimeDeriv: " << s.number_of_.TimeDeriv << std::endl;


		out << std::endl << " output locations:" << std::endl;
		out << "Functions " << s.output_locations_.Functions << std::endl;
		out << "Jacobian " << s.output_locations_.Jacobian << std::endl;

		if (s.HavePathVariable())
			out << "TimeDeriv " << s.output_locations_.TimeDeriv << std::endl;


		out << std::endl << " input locations:" << std::endl;
		out << "Variables " << s.input_locations_.Variables << std::endl;
		if (s.HavePathVariable())
			out << "Time " << s.input_locations_.Time << std::endl;



		out << std::endl << "true values of numbers: (number, location to downsample to)" << std::endl;
		for (auto const& x : s.true_values_of_numbers_)
		    out << *(x.first)  << ':' << x.second << std::endl;
		out << std::endl << std::endl;


		out << std::endl << "instructions: " << std::endl;
		for (size_t ii(0); ii<s.instructions_.size(); /*it's in the loop at access time*/){
			auto op = static_cast<Operation>(s.instructions_[ii++]);
			out << OpcodeToString(op) << "(";
			if (IsUnary(op))
				out << s.instructions_[ii++] << ") --> " << s.instructions_[ii++] << std::endl;

			else
				out << s.instructions_[ii++] << "," << s.instructions_[ii++] << ") --> " << s.instructions_[ii++] << std::endl;
		}




		auto& memory_dbl =  std::get<std::vector<dbl_complex>>(s.memory_);
		auto& memory_mpfr =  std::get<std::vector<mpfr_complex>>(s.memory_);

		out << "\nvariable values in dbl memory:\n";
		for (unsigned ii=0; ii<s.number_of_.Variables; ++ii){
			out << memory_dbl[s.input_locations_.Variables + ii] << " ";
		} 


		out << "\nvariable values in mpfr memory:\n";
		for (unsigned ii=0; ii<s.number_of_.Variables; ++ii){
			out << memory_mpfr[s.input_locations_.Variables + ii] << " ";
		} 

		if (s.HavePathVariable()){
			out << "\ntime value in dbl memory:\n";
				out << memory_dbl[s.input_locations_.Time] << " ";


			out << "\ntime value in mpfr memory:\n";
				out << memory_mpfr[s.input_locations_.Time] << " ";
		}

		out << std::endl << "full memory (double precision):" << std::endl;
		for (auto v: memory_dbl)
			out << v << ",";
		out << std::endl;

		
		out << std::endl << "full memory (mpfr precision):" << std::endl;
		for (auto v: memory_mpfr)
			out << v << ",";
		out << std::endl;



		return out;
	}


	template<typename NumT>
	void StraightLineProgram::Eval() const{

		auto& memory =  std::get<std::vector<NumT>>(memory_);


#ifndef BERTINI_DISABLE_PRECISION_CHECKS
		if (! std::is_same<NumT,dbl_complex>::value && Precision(memory[0])!=this->precision_){
			throw std::runtime_error("memory and SLP are out-of-sync WRT precision");
		}
#endif


		if (is_evaluated_)
			return;

		for (int ii = 0; ii<instructions_.size();/*the increment is done at end of loop depending on arity */) {
			//in the unary case the loop will increment by 3
			//binary: by 4

			switch (instructions_[ii]) {

				case Add:
					memory[this->instructions_[ii+3]] = memory[instructions_[ii+1]] + memory[instructions_[ii+2]];
					break;

				case Subtract:
					memory[this->instructions_[ii+3]] = memory[instructions_[ii+1]] - memory[instructions_[ii+2]];
					break;

				case Multiply:
					memory[this->instructions_[ii+3]] = memory[instructions_[ii+1]] * memory[instructions_[ii+2]];
					break;

				case Divide:
					memory[this->instructions_[ii+3]] = memory[instructions_[ii+1]] / memory[instructions_[ii+2]];
					break;

				case Power:
					memory[this->instructions_[ii+3]] = pow(memory[instructions_[ii+1]], memory[instructions_[ii+2]]);
					break;

				case IntPower:
					{
					memory[this->instructions_[ii+3]] = pow(memory[instructions_[ii+1]], this->integers_[instructions_[ii+2]]);
					break;
					}

				case Assign:
					memory[this->instructions_[ii+2]] = memory[instructions_[ii+1]];
					break;

				case Negate:
					memory[this->instructions_[ii+2]] = -(memory[instructions_[ii+1]]);
					break;

				case Sqrt:
					memory[this->instructions_[ii+2]] = sqrt(memory[instructions_[ii+1]]);
					break;

				case Log:
					memory[this->instructions_[ii+2]] = log(memory[instructions_[ii+1]]);
					break;

				case Exp:
					memory[this->instructions_[ii+2]] = exp(memory[instructions_[ii+1]]);
					break;

				case Sin:
					memory[this->instructions_[ii+2]] = sin(memory[instructions_[ii+1]]);
					break;

				case Cos:
					memory[this->instructions_[ii+2]] = cos(memory[instructions_[ii+1]]);
					break;

				case Tan:
					memory[this->instructions_[ii+2]] = tan(memory[instructions_[ii+1]]);
					break;

				case Asin:
					memory[this->instructions_[ii+2]] = asin(memory[instructions_[ii+1]]);
					break;

				case Acos:
					memory[this->instructions_[ii+2]] = acos(memory[instructions_[ii+1]]);
					break;

				case Atan:
					memory[this->instructions_[ii+2]] = atan(memory[instructions_[ii+1]]);
					break;

			} // switch for operation


			if (IsUnary(static_cast<Operation>(instructions_[ii]))) {
				ii = ii+3;

			}
			//in the binary case the loop will increment by 4
			else {
				ii = ii+4;
			}
		} // for loop around operations

		is_evaluated_ = true;
	}

	template void StraightLineProgram::Eval<dbl_complex>() const;
	template void StraightLineProgram::Eval<mpfr_complex>() const;


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


	//
	//
	//  implementer note:
	//
	// if you add another type to be visited, you must list it in TWO locations in the SLPCompiler type in the .hpp.   
	// 
	//



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

	void SLPCompiler::Visit(node::special_number::Pi const& n){
		this->DealWithNumber(n);
	}

	void SLPCompiler::Visit(node::special_number::E const& n){
		this->DealWithNumber(n);
	}


	void SLPCompiler::Visit(node::Jacobian const& n){
		throw std::runtime_error("unimplemented visit to node of type Jacobian");


	}

	void SLPCompiler::Visit(node::Differential const& n){
		throw std::runtime_error("unimplemented visit to node of type Differential");
	}





	void SLPCompiler::Visit(node::Function const & f){

		// put the location of the accepted node into memory, and copy into an output location.
		auto& n = f.entry_node();
		size_t location_entry;

		if (this->locations_encountered_nodes_.find(n) == this->locations_encountered_nodes_.end())
			n->Accept(*this); // think of calling Compile(n)

		location_entry = this->locations_encountered_nodes_[n];

		size_t location_this_node;
		if (this->locations_encountered_nodes_.find(f.shared_from_this()) == this->locations_encountered_nodes_.end()){
			location_this_node = next_available_complex_;
			locations_encountered_nodes_[f.shared_from_this()] = next_available_complex_++;
		}
		else
			location_this_node = locations_encountered_nodes_[f.shared_from_this()];


		slp_under_construction_.AddInstruction(Assign, location_entry, location_this_node);
	}


	// arithmetic
	void SLPCompiler::Visit(node::SumOperator const & n){

		// this loop
		// gets the locations of all the things we're going to add up.
		std::vector<size_t> operand_locations;
		for (auto& n : n.Operands()){

			if (this->locations_encountered_nodes_.find(n)==this->locations_encountered_nodes_.end())
				n->Accept(*this); // think of calling Compile(n)

			operand_locations.push_back(this->locations_encountered_nodes_[n]);
		}

		  
		const auto& signs = n.GetSigns();
		size_t prev_result_loc; // for tracking where the output of the previous iteration went

		// seed the loop.
		if (signs[0])
			prev_result_loc = operand_locations[0];
		else{
			slp_under_construction_.AddInstruction(Negate, operand_locations[0], next_available_complex_);
			prev_result_loc = next_available_complex_++;
		}
		

		// this loop
		// does the additions for the rest of the operands
		for (size_t ii{1}; ii<n.Operands().size(); ++ii){
			if (signs[ii])
				slp_under_construction_.AddInstruction(Add,prev_result_loc,operand_locations[ii],next_available_complex_);
			else
				slp_under_construction_.AddInstruction(Subtract,prev_result_loc,operand_locations[ii],next_available_complex_);

			prev_result_loc = next_available_complex_++;
		}

		this->locations_encountered_nodes_[n.shared_from_this()] =  prev_result_loc;

			// improved option?: we could do this in a way to minimize numerical error.  the obvious loop is not good for accumulation of error.  instead, use Pairwise summation
			// do pairs (0,1), (2,3), etc,
			// *then* add those temp vals together, until get to end.
			// see https://en.wikipedia.org/wiki/Pairwise_summation
	}







	void SLPCompiler::Visit(node::MultOperator const & n){
		


		// this loop
		// gets the locations of all the things we're going to add up.
		std::vector<size_t> operand_locations;
		for (auto& n : n.Operands()){

			if (this->locations_encountered_nodes_.find(n)==this->locations_encountered_nodes_.end())
				n->Accept(*this); // think of calling Compile(n)

			
			operand_locations.push_back(this->locations_encountered_nodes_[n]);
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
			auto location_one  = locations_encountered_nodes_[one];

			slp_under_construction_.AddInstruction(Divide, location_one, operand_locations[0], next_available_complex_);
			prev_result_loc = next_available_complex_++;
		}
		

		// this loop
		// does the additions for the rest of the operands
		for (size_t ii{1}; ii<n.Operands().size(); ++ii){
			if (mult_or_div[ii])
				slp_under_construction_.AddInstruction(Multiply,prev_result_loc,operand_locations[ii],next_available_complex_);
			else
				slp_under_construction_.AddInstruction(Divide,prev_result_loc,operand_locations[ii],next_available_complex_);

			prev_result_loc = next_available_complex_++;
		}

		this->locations_encountered_nodes_[n.shared_from_this()] =  prev_result_loc;

	}





	void SLPCompiler::Visit(node::IntegerPowerOperator const& n){
		
		IntT expo = n.exponent(); //integer

		// ensure we have the location of the base of the power operation.  it's a node at this point.
		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)
		auto location_operand = locations_encountered_nodes_[operand];



		if (this->locations_integers_.find(expo) == this->locations_integers_.end())
		{
			locations_integers_[expo] = slp_under_construction_.integers_.size();
			slp_under_construction_.integers_.push_back(expo);
		}

		auto location_exponent  = locations_integers_[expo]; // this is a map lookup

		
		this->locations_encountered_nodes_[n.shared_from_this()] = next_available_complex_;
		slp_under_construction_.AddInstruction(IntPower,location_operand,location_exponent, next_available_complex_++);
	}


	void SLPCompiler::Visit(node::PowerOperator const& n){
		
		//get location of base and power then add instruction

		const auto& base = n.GetBase();
		const auto& exponent = n.GetExponent();

		if (this->locations_encountered_nodes_.find(base) == this->locations_encountered_nodes_.end())
			base->Accept(*this); // think of calling Compile(n)

		if (this->locations_encountered_nodes_.find(exponent) == this->locations_encountered_nodes_.end())
			exponent->Accept(*this); // think of calling Compile(n)

		auto loc_base = locations_encountered_nodes_[base];
		auto loc_exponent = locations_encountered_nodes_[exponent];

		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Power, loc_base, loc_exponent, next_available_complex_++);



	}

	void SLPCompiler::Visit(node::ExpOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Exp,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::LogOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Log,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::NegateOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Negate,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::SqrtOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Sqrt,location_operand, next_available_complex_++);
	}


	// the trig operators
	void SLPCompiler::Visit(node::SinOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Sin,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::ArcSinOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Asin,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::CosOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Cos,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::ArcCosOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Acos,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::TanOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Tan,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::ArcTanOperator const& n){
		

		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this); // think of calling Compile(n)

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[n.shared_from_this()] =  next_available_complex_;
		slp_under_construction_.AddInstruction(Atan,location_operand, next_available_complex_++);
	}



	SLP SLPCompiler::Compile(System const& sys){
		this->Clear();

		


		this->slp_under_construction_.precision_ = DefaultPrecision();

		// deal with variables
		

			// 1. ADD VARIABLES

		slp_under_construction_.input_locations_.Variables = next_available_mem_;

		auto variable_ordering = sys.VariableOrdering();
		for (auto v: variable_ordering){
			locations_encountered_nodes_[ v ] = next_available_complex_++;
		}
		slp_under_construction_.number_of_.Variables = variable_ordering.size();
		// slp_under_construction_.input_locations_.Variables = variable_counter;

			// deal with path variable
		if (sys.HavePathVariable())
		{
				// do this action only if the system has a path variable defined
			slp_under_construction_.input_locations_.Time = next_available_complex_;
			locations_encountered_nodes_[ sys.GetPathVariable() ] = next_available_complex_++;
			slp_under_construction_.has_path_variable_ = true;
		}


		
			// make space for natural functions and derivatives.  we omit the patches.
			// 3. ADD FUNCTIONS
		slp_under_construction_.number_of_.Functions = sys.NumNaturalFunctions();
		slp_under_construction_.output_locations_.Functions = next_available_complex_;
		for (auto f: sys.GetNaturalFunctions())
			locations_encountered_nodes_[f] = next_available_complex_++;





		
		// always have space derivatives

		auto ds_dx = sys.GetSpaceDerivatives(); // a linear object, so can just run down the object
		slp_under_construction_.number_of_.Jacobian = ds_dx.size();
		slp_under_construction_.output_locations_.Jacobian = next_available_complex_;
		for (auto n: ds_dx)
			locations_encountered_nodes_[n] = next_available_complex_++;


		// sometimes have time derivatives
		if (sys.HavePathVariable()) {
			
			auto ds_dt = sys.GetTimeDerivatives();  // a linear object, so can just run down the object
			slp_under_construction_.number_of_.TimeDeriv = ds_dt.size();
			slp_under_construction_.output_locations_.TimeDeriv = next_available_complex_;
			for (auto n: ds_dt)
				locations_encountered_nodes_[n] = next_available_complex_++;
		}






		
		for (auto f: sys.GetNaturalFunctions())
		{
			f->Accept(*this);

			// post visit function
			/* code */
		}



		
		// always do derivatives with respect to space variables
		// 4. ADD SPACE VARIABLE DERIVATIVES
		for (auto n: ds_dx)
			n->Accept(*this);




		// sometimes have time derivatives
		if (sys.HavePathVariable()) {
			
			// we need derivatives with respect to time only if the system has a path variable defined
			// 5. ADD TIME VARIABLE DERIVATIVES

			auto ds_dt = sys.GetTimeDerivatives();  // a linear object, so can just run down the object
			for (auto n: ds_dt)
				n->Accept(*this);
		}

		slp_under_construction_.GetMemory<dbl_complex>().resize(next_available_complex_);
		slp_under_construction_.GetMemory<mpfr_complex>().resize(next_available_complex_);

		slp_under_construction_.CopyNumbersIntoMemory<dbl_complex>();
		slp_under_construction_.CopyNumbersIntoMemory<mpfr_complex>();

		return slp_under_construction_;
	}

	void SLPCompiler::Clear(){
		next_available_complex_ = 0;
		next_available_int_ = 0;

		locations_encountered_nodes_.clear();
		slp_under_construction_ = SLP();
	}
}
