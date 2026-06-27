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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire
// michael mumm, university of wisconsin eau claire

#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/system/blocks/polynomial_block.hpp"

#include <boost/math/constants/constants.hpp>



BOOST_CLASS_EXPORT(bertini::StraightLineProgram);

// SLP Stuff
namespace bertini{
	using namespace bertini::node;

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
		throw std::runtime_error("unrecognized operation in OpcodeToString");
	}


	// Produce a constant's value directly from its exact recipe --- no function-tree node, no node
	// evaluation (ADR-0027).  These mirror the number nodes' FreshEval_d / FreshEval_mp exactly
	// (Integer/Complex/Rational in number.cpp; Pi/E in special_number.cpp) so the compiled program
	// evaluates bit-for-bit identically to the old node-backed path, at the ambient working
	// precision (ThreadPrecision).
	template<>
	complex_dbl ConstantRecipe::Produce<complex_dbl>() const {
		switch (kind) {
			case Kind::Integer:  return complex_dbl(double(int_value), 0);
			case Kind::Rational: return complex_dbl(double(rat_real), double(rat_imag));
			case Kind::Complex:    return complex_dbl(float_value);
			case Kind::Pi:       return complex_dbl(boost::math::constants::pi<double>(), 0);
			case Kind::E:        return complex_dbl(exp(1.0), 0.0);
		}
		throw std::runtime_error("unrecognized ConstantRecipe kind in Produce<complex_dbl>");
	}

	template<>
	complex_mp ConstantRecipe::Produce<complex_mp>() const {
		using boost::multiprecision::mpfr_float;
		switch (kind) {
			case Kind::Integer:  return complex_mp(int_value, 0, ThreadPrecision());
			case Kind::Rational: return complex_mp(real_mp(rat_real, ThreadPrecision()), real_mp(rat_imag, ThreadPrecision()));
			case Kind::Complex:    return complex_mp(float_value, ThreadPrecision());
			case Kind::Pi:       return complex_mp(boost::math::constants::pi<real_mp>());
			case Kind::E:        return complex_mp(real_mp(exp(real_mp(1))));
		}
		throw std::runtime_error("unrecognized ConstantRecipe kind in Produce<complex_mp>");
	}

	// Real companion of Produce (ADR-0034): the real value of a constant whose slot was inferred
	// NumType::Real.  Only the real part is used; IsReal() guarantees the imaginary part is zero, so
	// this matches the real part of Produce<complex>() bit-for-bit.
	template<>
	real_dbl ConstantRecipe::ProduceReal<real_dbl>() const {
		switch (kind) {
			case Kind::Integer:  return double(int_value);
			case Kind::Rational: return double(rat_real);
			case Kind::Complex:  return double(float_value.real());
			case Kind::Pi:       return boost::math::constants::pi<double>();
			case Kind::E:        return exp(1.0);
		}
		throw std::runtime_error("unrecognized ConstantRecipe kind in ProduceReal<real_dbl>");
	}

	template<>
	real_mp ConstantRecipe::ProduceReal<real_mp>() const {
		switch (kind) {
			case Kind::Integer:  return real_mp(int_value, ThreadPrecision());
			case Kind::Rational: return real_mp(rat_real, ThreadPrecision());
			case Kind::Complex:  return real_mp(float_value.real(), ThreadPrecision());
			case Kind::Pi:       return boost::math::constants::pi<real_mp>();
			case Kind::E:        return real_mp(exp(real_mp(1)));
		}
		throw std::runtime_error("unrecognized ConstantRecipe kind in ProduceReal<real_mp>");
	}


	// the constructor
	StraightLineProgram::StraightLineProgram(System const& sys){
		SLPCompiler compiler;

		*this = compiler.Compile(sys);
	}

	void StraightLineProgram::precision(unsigned new_precision) const{

		if (new_precision==memory_.precision_){
			return;
		}
		else{
			auto& cmem = memory_.Get<complex_mp>();
			auto& rmem = memory_.Get<real_mp>();

			// Refill the constants from their exact recipes (no node evaluation) into the bank their
			// NumType selects, then normalize every slot of both banks to the new precision.
			for (auto const& c : program_->constant_recipes_){
				if (program_->slot_numtype_[c.slot] == NumType::Real)
					rmem[c.slot] = c.ProduceReal<real_mp>();
				else
					cmem[c.slot] = c.Produce<complex_mp>();
			}

			for (auto& n : cmem)
				Precision(n, new_precision);
			for (auto& n : rmem)
				Precision(n, new_precision);

			memory_.precision_ = new_precision;
		}


	}



	template<typename NumT>
	void StraightLineProgram::CopyNumbersIntoMemory() const
	{
		using RealT = typename NumTraits<NumT>::Real;
		constexpr bool is_mp = std::is_same<NumT,complex_mp>::value;
		for (auto const& c : program_->constant_recipes_){
			if (program_->slot_numtype_[c.slot] == NumType::Real){
				memory_.Get<RealT>()[c.slot] = c.ProduceReal<RealT>();
				if constexpr (is_mp)
					Precision(memory_.Get<RealT>()[c.slot], memory_.precision_);
			}
			else{
				memory_.Get<NumT>()[c.slot] = c.Produce<NumT>();
				if constexpr (is_mp)
					Precision(memory_.Get<NumT>()[c.slot], memory_.precision_);
			}
		}
	}

	template void StraightLineProgram::CopyNumbersIntoMemory<complex_dbl>() const;
	template void StraightLineProgram::CopyNumbersIntoMemory<complex_mp>() const;


	void StraightLineProgram::SetupMemory()
	{
		// Re-seed the thread-local default precision from the program's own precision before
		// growing the complex_mp memory block. resize() default-constructs each new element via
		// mpfr_init2(x, thread_default_precision()); on Boost 1.87 that value can be 0 on fresh
		// threads, which aborts. DefaultPrecision() sets both the static and thread-local defaults
		// so the default-constructed slots are valid.
		DefaultPrecision(memory_.precision_);

		// adjust the sizes of the memory blocks to match the number expected via compilation.
		// All four banks (real/complex x dbl/mp) are sized to the slot count; a slot lives in exactly
		// one bank per its NumType, so the off-type banks hold default, never-read entries for now.
		memory_.Get<real_dbl>().resize(program_->num_slots_);
		memory_.Get<complex_dbl>().resize(program_->num_slots_);
		memory_.Get<real_mp>().resize(program_->num_slots_);
		memory_.Get<complex_mp>().resize(program_->num_slots_);

		// downsample to get ready for evaluation
		CopyNumbersIntoMemory<complex_dbl>();
		CopyNumbersIntoMemory<complex_mp>();
	}




	void SLPProgram::AddInstruction(Operation binary_op, size_t in_loc1, size_t in_loc2, size_t out_loc){
		this->instructions_.push_back(binary_op);
		this->instructions_.push_back(in_loc1);
		this->instructions_.push_back(in_loc2);
		this->instructions_.push_back(out_loc);
	}



	void SLPProgram::AddInstruction(Operation unary_op, size_t in_loc, size_t out_loc){
		this->instructions_.push_back(unary_op);
		this->instructions_.push_back(in_loc);
		this->instructions_.push_back(out_loc);

	}

	void SLPProgram::AddConstant(ConstantRecipe recipe){
		this->constant_recipes_.push_back(std::move(recipe));
	}



	std::ostream& operator <<(std::ostream& out, const StraightLineProgram & s){
		auto const& prog = *s.program_;

		out << "\n\n#fns: " << s.NumFunctions() << " #vars: " << s.NumVariables() << std::endl;
		out << "have path variable: " << s.HavePathVariable() << std::endl;

		out << std::endl << "numbers of things:" << std::endl;
		out << "Functions: " << prog.number_of_.Functions << std::endl;
		out << "Variables: " << prog.number_of_.Variables << std::endl;
		out << "Jacobian: " << prog.number_of_.Jacobian << std::endl;

		if (s.HavePathVariable())
			out << "TimeDeriv: " << prog.number_of_.TimeDeriv << std::endl;


		out << std::endl << " output locations:" << std::endl;
		out << "Functions " << prog.output_locations_.Functions << std::endl;
		out << "Jacobian " << prog.output_locations_.Jacobian << std::endl;

		if (s.HavePathVariable())
			out << "TimeDeriv " << prog.output_locations_.TimeDeriv << std::endl;


		out << std::endl << " input locations:" << std::endl;
		out << "Variables " << prog.input_locations_.Variables << std::endl;
		if (s.HavePathVariable())
			out << "Time " << prog.input_locations_.Time << std::endl;



		out << std::endl << "constants: (kind, location to downsample to)" << std::endl;
		for (auto const& c : prog.constant_recipes_)
		    out << static_cast<int>(c.kind)  << ':' << c.slot << std::endl;
		out << std::endl << std::endl;


		out << std::endl << "instructions: " << std::endl;
		for (size_t ii(0); ii<prog.instructions_.size(); /*it's in the loop at access time*/){
			auto op = static_cast<Operation>(prog.instructions_[ii++]);
			out << OpcodeToString(op) << "(";
			if (IsUnary(op)){
				auto operand = prog.instructions_[ii++];
				auto result = prog.instructions_[ii++];
				out << operand << ") --> " << result << std::endl;
			}
			else{
				auto operand1 = prog.instructions_[ii++];
				auto operand2 = prog.instructions_[ii++];
				auto result = prog.instructions_[ii++];
				out << operand1 << "," << operand2 << ") --> " << result << std::endl;
			}
		}




		auto& memory_dbl =  s.memory_.Get<complex_dbl>();
		auto& memory_mpfr =  s.memory_.Get<complex_mp>();

		out << "\nvariable values in complex_dbl memory:\n";
		for (unsigned ii=0; ii<prog.number_of_.Variables; ++ii){
			out << memory_dbl[prog.input_locations_.Variables + ii] << " ";
		}


		out << "\nvariable values in mpfr memory:\n";
		for (unsigned ii=0; ii<prog.number_of_.Variables; ++ii){
			out << memory_mpfr[prog.input_locations_.Variables + ii] << " ";
		}

		if (s.HavePathVariable()){
			out << "\ntime value in complex_dbl memory:\n";
				out << memory_dbl[prog.input_locations_.Time] << " ";


			out << "\ntime value in mpfr memory:\n";
				out << memory_mpfr[prog.input_locations_.Time] << " ";
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
	void SLPProgram::Eval(SLPMemory& mem) const{

		// Two banks at this precision: the complex bank (NumT) and its real companion (ADR-0034).
		// A slot's value lives in exactly one, selected by slot_numtype_.  An all-Complex program
		// only ever touches `cplx`, so this is identical to the pre-tier path.
		using RealT = typename NumTraits<NumT>::Real;
		auto& cplx = mem.Get<NumT>();
		auto& real = mem.Get<RealT>();

		// Bring the std math functions into scope so the REAL (double) path resolves them: a plain
		// double has no associated namespace, so unqualified sin/pow/... would otherwise bind to the
		// function-tree node operators.  ADL still finds Boost's overloads for real_mp / complex_mp
		// and std's for complex_dbl, so all four numeric types resolve correctly.
		using std::pow;  using std::sqrt; using std::log;  using std::exp;
		using std::sin;  using std::cos;  using std::tan;
		using std::asin; using std::acos; using std::atan;

		// The tier dispatch constructs complex temporaries during eval (promoting a real operand to
		// complex in Power/sqrt/log/...).  Boost inits such a new mpc at the *thread default* precision,
		// which is not guaranteed to be the working precision on this path -- if it is 0, mpc_init2
		// aborts.  Pin the thread-local default to the working precision (thread-safe; no global write).
		// (The old eval never constructed mp temporaries, so it didn't need this.)
		if constexpr (!std::is_same<NumT,complex_dbl>::value)
			SetThreadPrecision(mem.precision_);

		// Re-tag a freshly written complex slot to the working precision.  Boost.Multiprecision has a
		// bug in mixed mpfr_float/mpc_complex expression templates: `real_mp / complex_mp` computes the
		// correct value but tags the result precision 0 (multiplication is unaffected), which later
		// aborts when that slot is copied (mpc_init2 with precision 0).  Minimal bertini-free repro:
		// two operands at precision 40, `out = r / z` gives out.precision()==0.  Re-tagging restores it.
		// No-op for complex_dbl (std::complex has no precision tag), so the hot double path keeps native
		// mixed arithmetic with zero overhead.
		auto retag = [&](size_t o){
			if constexpr (!std::is_same<NumT,complex_dbl>::value) Precision(cplx[o], mem.precision_);
		};


#ifndef BERTINI_DISABLE_PRECISION_CHECKS
		if (! std::is_same<NumT,complex_dbl>::value && Precision(cplx[0])!=mem.precision_){
			throw std::runtime_error("memory and SLP are out-of-sync WRT precision");
		}
#endif


		if (mem.is_evaluated_)
			return;

		// If the frozen prologue's results are already valid for this number type (double constants
		// never change; mpfr constants are valid while the working precision is unchanged), skip it
		// and re-run only the live segment, reusing the frozen slots already in memory.
		bool frozen_valid;
		if constexpr (std::is_same<NumT,complex_dbl>::value)
			frozen_valid = mem.frozen_valid_dbl_;
		else
			frozen_valid = (mem.frozen_valid_mp_precision_ == mem.precision_);

		const size_t loop_start = frozen_valid ? first_live_instruction_ : 0;

		auto is_real = [&](size_t s){ return slot_numtype_[s] == NumType::Real; };

		// Binary arithmetic: read each operand from its NumType's bank and write the result to the
		// result slot's bank.  When the result is Real, inference guarantees both operands are Real
		// (R+R, R-R, R*R, R/R), so we use the cheap real path; otherwise we use native mixed
		// arithmetic (real * complex is ~half the work of complex * complex), promoting nothing.
		// `is_div` marks Divide: only `real_mp / complex_mp` trips the Boost precision bug (and only
		// when the complex value's imaginary part is 0, which is data-dependent), so we re-tag exactly
		// that one branch -- the common mixed multiply/add/subtract need no fix-up.
		auto binop = [&](size_t i1, size_t i2, size_t o, bool is_div, auto fn){
			if (is_real(o))                          real[o] = fn(real[i1], real[i2]);  // all-real, cheap
			else if (!is_real(i1) && !is_real(i2))   cplx[o] = fn(cplx[i1], cplx[i2]);  // pure complex
			else if (is_real(i1)) {                                                     // real (op) complex
				cplx[o] = fn(real[i1], cplx[i2]);
				if (is_div) retag(o);
			}
			else                                     cplx[o] = fn(cplx[i1], real[i2]);  // complex (op) real
		};

		// Type-preserving unary (result NumType == operand NumType): negate/copy/exp/sin/cos/tan/atan.
		auto un_preserve = [&](size_t i, size_t o, auto fn){
			if (is_real(o)) real[o] = fn(real[i]); else cplx[o] = fn(cplx[i]);
		};

		// Escape unary (result is Complex; a real operand is promoted so the complex branch is taken):
		// sqrt/log/asin/acos, which can leave ℝ for some real inputs.
		auto un_escape = [&](size_t i, size_t o, auto fn){
			// a real operand is promoted by a single converting construction (at the working precision),
			// not a mixed expression template, so the result is tagged correctly -- no re-tag needed.
			cplx[o] = is_real(i) ? fn(NumT(real[i])) : fn(cplx[i]);
		};

		for (size_t ii = loop_start; ii<instructions_.size();/*the increment is done at end of loop depending on arity */) {
			//in the unary case the loop will increment by 3
			//binary: by 4

			const size_t a = instructions_[ii+1], b = instructions_[ii+2], c = instructions_[ii+3];

			switch (instructions_[ii]) {

				case Add:      binop(a, b, c, false, [](auto const& x, auto const& y){ return x + y; }); break;
				case Subtract: binop(a, b, c, false, [](auto const& x, auto const& y){ return x - y; }); break;
				case Multiply: binop(a, b, c, false, [](auto const& x, auto const& y){ return x * y; }); break;
				case Divide:   binop(a, b, c, true,  [](auto const& x, auto const& y){ return x / y; }); break;

				case Power: {
					// general a^b (slot exponent); result is Complex, so promote any real operand to
					// complex (a single converting construction at the working precision, not a mixed
					// expression template), then a pure-complex pow -- so the result is tagged correctly.
					const NumT base = is_real(a) ? NumT(real[a]) : cplx[a];
					const NumT expo = is_real(b) ? NumT(real[b]) : cplx[b];
					cplx[c] = pow(base, expo);
					break;
				}

				case IntPower:
					// in2 (b) is an index into integers_, not a slot.  real^int -> real, complex^int -> complex.
					if (is_real(c)) real[c] = pow(real[a], this->integers_[b]);
					else            cplx[c] = pow(cplx[a], this->integers_[b]);
					break;

				case Assign:   un_preserve(a, b, [](auto const& x){ return x; }); break;
				case Negate:   un_preserve(a, b, [](auto const& x){ return -x; }); break;
				case Exp:      un_preserve(a, b, [](auto const& x){ return exp(x); }); break;
				case Sin:      un_preserve(a, b, [](auto const& x){ return sin(x); }); break;
				case Cos:      un_preserve(a, b, [](auto const& x){ return cos(x); }); break;
				case Tan:      un_preserve(a, b, [](auto const& x){ return tan(x); }); break;
				case Atan:     un_preserve(a, b, [](auto const& x){ return atan(x); }); break;

				case Sqrt:     un_escape(a, b, [](auto const& x){ return sqrt(x); }); break;
				case Log:      un_escape(a, b, [](auto const& x){ return log(x); }); break;
				case Asin:     un_escape(a, b, [](auto const& x){ return asin(x); }); break;
				case Acos:     un_escape(a, b, [](auto const& x){ return acos(x); }); break;

			} // switch for operation


			if (IsUnary(static_cast<Operation>(instructions_[ii]))) {
				ii = ii+3;

			}
			//in the binary case the loop will increment by 4
			else {
				ii = ii+4;
			}
		} // for loop around operations

		// A full run (from instruction 0) has just refreshed the frozen prologue at this precision.
		if (!frozen_valid)
		{
			if constexpr (std::is_same<NumT,complex_dbl>::value)
				mem.frozen_valid_dbl_ = true;
			else
				mem.frozen_valid_mp_precision_ = mem.precision_;
		}

		mem.is_evaluated_ = true;
	}

	template void SLPProgram::Eval<complex_dbl>(SLPMemory&) const;
	template void SLPProgram::Eval<complex_mp>(SLPMemory&) const;


	namespace {
		// The NumType of a binary op's result given its operands' NumTypes (ADR-0034).  Real only when
		// the result is GUARANTEED real for all inputs; otherwise Complex (the safe escape).  IntPower
		// is handled separately (its exponent is an integer index, and result = base NumType).
		NumType BinaryResultNumType(Operation op, NumType a, NumType b)
		{
			const bool both_real = (a == NumType::Real && b == NumType::Real);
			switch (op)
			{
				// R+R, R-R, R*R are real; R/R is real (real/real stays in ℝ); any complex operand -> complex.
				case Add: case Subtract: case Multiply: case Divide:
					return both_real ? NumType::Real : NumType::Complex;
				// general a^b (slot exponent): negative real base to a non-integer power leaves ℝ.
				case Power:
				default:
					return NumType::Complex;
			}
		}

		// The NumType of a unary op's result given its operand's NumType.
		NumType UnaryResultNumType(Operation op, NumType a)
		{
			switch (op)
			{
				case Negate: case Assign:                          return a;          // sign/copy preserve
				case Exp: case Sin: case Cos: case Tan: case Atan: return a;          // real-valued for real input
				case Sqrt: case Log: case Asin: case Acos:         return NumType::Complex;  // can leave ℝ
				default:                                           return NumType::Complex;
			}
		}
	} // anonymous namespace

	bool SLPProgram::tiers_enabled_ = true;

	void SLPProgram::ComputeSlotNumTypes()
	{
		// Seed: constant slots from their recipe's real-ness; everything else (variables, the path
		// variable, and as-yet-unwritten temporaries) starts Complex.  Then one forward pass over the
		// dependency-ordered tape propagates the NumType of each instruction's result (same shape as
		// PartitionInstructions' frozenness pass).
		slot_numtype_.assign(num_slots_, NumType::Complex);
		if (!tiers_enabled_)  // A/B baseline: force the pre-tier all-complex evaluation
			return;
		for (auto const& c : constant_recipes_)
			slot_numtype_[c.slot] = c.IsReal() ? NumType::Real : NumType::Complex;

		for (size_t ii = 0; ii < instructions_.size(); )
		{
			const auto op = static_cast<Operation>(instructions_[ii]);
			if (IsUnary(op))
			{
				slot_numtype_[instructions_[ii + 2]] =
					UnaryResultNumType(op, slot_numtype_[instructions_[ii + 1]]);
				ii += 3;
			}
			else if (op == IntPower)
			{
				// in2 is an index into integers_, not a slot; real^int = real, complex^int = complex.
				slot_numtype_[instructions_[ii + 3]] = slot_numtype_[instructions_[ii + 1]];
				ii += 4;
			}
			else
			{
				slot_numtype_[instructions_[ii + 3]] =
					BinaryResultNumType(op, slot_numtype_[instructions_[ii + 1]], slot_numtype_[instructions_[ii + 2]]);
				ii += 4;
			}
		}
	}


	void SLPProgram::PartitionInstructions()
	{
		// A memory slot is "frozen" if its value depends only on frozen inputs.  Seed: the literal
		// numbers (Integer/Complex/Rational and Pi/E, all in true_values_of_numbers_) are frozen; the
		// variable and time slots are live.  Then a single forward pass propagates frozenness: an
		// instruction is frozen iff all its input slots are frozen, and it freezes its output slot.
		const size_t num_slots = num_slots_;
		std::vector<bool> slot_frozen(num_slots, false);
		for (auto const& c : constant_recipes_)
			slot_frozen[c.slot] = true;

		struct Instr { size_t off; size_t len; bool frozen; };
		std::vector<Instr> parsed;

		for (size_t ii = 0; ii < instructions_.size(); )
		{
			const auto op = static_cast<Operation>(instructions_[ii]);
			const bool unary = IsUnary(op);
			const size_t len = unary ? 3 : 4;

			bool frozen;
			size_t out;
			if (unary)
			{
				out = instructions_[ii + 2];
				frozen = slot_frozen[instructions_[ii + 1]];
			}
			else if (op == IntPower)
			{
				// The second operand of IntPower is an index into integers_, not a memory slot; the
				// exponent is a literal, so frozenness depends only on the base slot.
				out = instructions_[ii + 3];
				frozen = slot_frozen[instructions_[ii + 1]];
			}
			else
			{
				out = instructions_[ii + 3];
				frozen = slot_frozen[instructions_[ii + 1]] && slot_frozen[instructions_[ii + 2]];
			}

			slot_frozen[out] = frozen;
			parsed.push_back({ii, len, frozen});
			ii += len;
		}

		// Stable partition: frozen instructions first (preserving relative order), then live ones.
		// This is dependency-safe because no live instruction is an input to a frozen one.
		std::vector<size_t> reordered;
		reordered.reserve(instructions_.size());
		size_t frozen_words = 0;
		for (auto const& I : parsed)
			if (I.frozen)
			{
				reordered.insert(reordered.end(), instructions_.begin() + static_cast<std::ptrdiff_t>(I.off), instructions_.begin() + static_cast<std::ptrdiff_t>(I.off) + static_cast<std::ptrdiff_t>(I.len));
				frozen_words += I.len;
			}
		for (auto const& I : parsed)
			if (!I.frozen)
				reordered.insert(reordered.end(), instructions_.begin() + static_cast<std::ptrdiff_t>(I.off), instructions_.begin() + static_cast<std::ptrdiff_t>(I.off) + static_cast<std::ptrdiff_t>(I.len));

		instructions_ = std::move(reordered);
		first_live_instruction_ = frozen_words;
	}

}







// stuff for SLPCompiler
namespace bertini{
	using SLP = StraightLineProgram;


	// Build an exact ConstantRecipe straight from a number node's true value --- no node evaluation
	// (ADR-0027).  One overload per concrete constant kind.
	namespace {
		ConstantRecipe RecipeFor(node::Integer const& n){
			ConstantRecipe r; r.kind = ConstantRecipe::Kind::Integer; r.int_value = n.GetValue(); return r;
		}
		ConstantRecipe RecipeFor(node::Rational const& n){
			ConstantRecipe r; r.kind = ConstantRecipe::Kind::Rational;
			r.rat_real = n.GetValueReal(); r.rat_imag = n.GetValueImag(); return r;
		}
		ConstantRecipe RecipeFor(node::Complex const& n){
			ConstantRecipe r; r.kind = ConstantRecipe::Kind::Complex; r.float_value = n.GetValue(); return r;
		}
		ConstantRecipe RecipeFor(node::special_number::Pi const&){
			ConstantRecipe r; r.kind = ConstantRecipe::Kind::Pi; return r;
		}
		ConstantRecipe RecipeFor(node::special_number::E const&){
			ConstantRecipe r; r.kind = ConstantRecipe::Kind::E; return r;
		}
	}


	void SLPCompiler::RegisterConstant(Nd const& nd, ConstantRecipe recipe){
		recipe.slot = next_available_complex_;
		program_under_construction_.AddConstant(std::move(recipe));
		locations_encountered_nodes_[nd] = next_available_complex_++;
	}


	void SLPCompiler::Visit(node::Variable const& n){
		// A system's variables are all pre-registered before its function trees are compiled, so
		// reaching this Visit means a function references a variable that is not in the system's
		// variable ordering.  That is unsupported: to bake a constant into a function, build it with
		// a literal (Complex / Integer / Rational), not a variable removed from the ordering.
		throw std::runtime_error("SLP compile: a function references the variable '" + n.name() +
			"', which is not in the system's variable ordering");
	}


	//
	//
	//  implementer note:
	//
	// if you add another type to be visited, you must list it in TWO locations in the SLPCompiler type in the .hpp.
	//
	//



	void SLPCompiler::Visit(node::Integer const& n){
		this->RegisterConstant(n.shared_from_this(), RecipeFor(n));
	}

	void SLPCompiler::Visit(node::Complex const& n){
		this->RegisterConstant(n.shared_from_this(), RecipeFor(n));
	}

	void SLPCompiler::Visit(node::Rational const& n){
		this->RegisterConstant(n.shared_from_this(), RecipeFor(n));
	}

	void SLPCompiler::Visit(node::special_number::Pi const& n){
		this->RegisterConstant(n.shared_from_this(), RecipeFor(n));
	}

	void SLPCompiler::Visit(node::special_number::E const& n){
		this->RegisterConstant(n.shared_from_this(), RecipeFor(n));
	}


	void SLPCompiler::Visit(node::Differential const& /*n*/){
		throw std::runtime_error("unimplemented visit to node of type Differential");
	}





	void SLPCompiler::Visit(node::NamedExpression const & f){
		// A named subexpression appearing inside a tree (a = x^2+y^2, used elsewhere): compute the
		// entry once and copy its value into the NamedExpression's own slot, so every reference to
		// the name shares that one result.  (Same wiring as an embedded Function.)
		const std::shared_ptr<node::Node>& n = f.EntryNode();
		const std::shared_ptr<const node::NamedExpression> f_as_ptr = std::dynamic_pointer_cast<node::NamedExpression const>(f.shared_from_this());

		if (this->locations_encountered_nodes_.find(n) == this->locations_encountered_nodes_.end())
			n->Accept(*this);
		size_t location_entry = this->locations_encountered_nodes_[n];

		size_t location_this_node;
		if (this->locations_encountered_nodes_.find(f_as_ptr) == this->locations_encountered_nodes_.end()){
			location_this_node = next_available_complex_;
			locations_encountered_nodes_[f_as_ptr] = next_available_complex_++;
		}
		else
			location_this_node = locations_encountered_nodes_[f_as_ptr];

		program_under_construction_.AddInstruction(Assign, location_entry, location_this_node);
	}


	// arithmetic
	void SLPCompiler::Visit(node::SumOperator const & n){
		const std::shared_ptr<const node::SumOperator> as_ptr = std::dynamic_pointer_cast<node::SumOperator const>(n.shared_from_this());

		// this loop
		// gets the locations of all the things we're going to add up.
		std::vector<size_t> operand_locations;
		for (auto& n : n.Operands()){

			if (this->locations_encountered_nodes_.find(n)==this->locations_encountered_nodes_.end())
				n->Accept(*this);

			operand_locations.push_back(this->locations_encountered_nodes_[n]);
		}


		const auto& signs = n.GetSigns();
		size_t prev_result_loc; // for tracking where the output of the previous iteration went

		// seed the loop.
		if (signs[0])
			prev_result_loc = operand_locations[0];
		else{
			program_under_construction_.AddInstruction(Negate, operand_locations[0], next_available_complex_);
			prev_result_loc = next_available_complex_++;
		}


		// this loop
		// does the additions for the rest of the operands
		for (size_t ii{1}; ii<n.Operands().size(); ++ii){
			if (signs[ii])
				program_under_construction_.AddInstruction(Add,prev_result_loc,operand_locations[ii],next_available_complex_);
			else
				program_under_construction_.AddInstruction(Subtract,prev_result_loc,operand_locations[ii],next_available_complex_);

			prev_result_loc = next_available_complex_++;
		}

		this->locations_encountered_nodes_[as_ptr] =  prev_result_loc;

			// improved option?: we could do this in a way to minimize numerical error.  the obvious loop is not good for accumulation of error.  instead, use Pairwise summation
			// do pairs (0,1), (2,3), etc,
			// *then* add those temp vals together, until get to end.
			// see https://en.wikipedia.org/wiki/Pairwise_summation
	}







	void SLPCompiler::Visit(node::MultOperator const & n){

		const std::shared_ptr<const node::MultOperator> as_ptr = std::dynamic_pointer_cast<node::MultOperator const>(n.shared_from_this());

		// this loop
		// gets the locations of all the things we're going to add up.
		std::vector<size_t> operand_locations;
		for (auto& operand : n.Operands()){

			if (this->locations_encountered_nodes_.find(operand)==this->locations_encountered_nodes_.end())
			{
				operand->Accept(*this);
			}

			operand_locations.push_back(this->locations_encountered_nodes_[operand]);
		}


		const auto& mult_or_div = n.GetMultOrDiv();// true is multiply and false is divide

		size_t prev_result_loc; // for tracking where the output of the previous iteration went

		// seed the loop.
		if (mult_or_div[0])
			prev_result_loc = operand_locations[0];
		else{
			// this case is reciprocation of the first operand

			// this code sucks.  really, there should be a bank of integers that we pull from, instead of many copies of the same integer.
			auto one = Integer::Make(1);
			this->RegisterConstant(one, RecipeFor(*one));
			auto location_one  = locations_encountered_nodes_[one];

			program_under_construction_.AddInstruction(Divide, location_one, operand_locations[0], next_available_complex_);
			prev_result_loc = next_available_complex_++;
		}


		// this loop
		// does the additions for the rest of the operands
		for (size_t ii{1}; ii<n.Operands().size(); ++ii){
			if (mult_or_div[ii])
				program_under_construction_.AddInstruction(Multiply,prev_result_loc,operand_locations[ii],next_available_complex_);
			else
				program_under_construction_.AddInstruction(Divide,prev_result_loc,operand_locations[ii],next_available_complex_);

			prev_result_loc = next_available_complex_++;
		}

		this->locations_encountered_nodes_[as_ptr] =  prev_result_loc;

	}





	void SLPCompiler::Visit(node::IntegerPowerOperator const& n){
		auto as_ptr = std::dynamic_pointer_cast<node::IntegerPowerOperator const>(n.shared_from_this());

		IntT expo = n.exponent(); //integer

		// ensure we have the location of the base of the power operation.  it's a node at this point.
		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
		{
			operand->Accept(*this);
		}

		auto location_operand = locations_encountered_nodes_[operand];



		if (this->locations_integers_.find(expo) == this->locations_integers_.end())
		{
			locations_integers_[expo] = program_under_construction_.integers_.size();
			program_under_construction_.integers_.push_back(expo);
		}


		auto location_exponent  = locations_integers_[expo]; // this is a map lookup


		this->locations_encountered_nodes_[as_ptr] = next_available_complex_;
		program_under_construction_.AddInstruction(IntPower,location_operand,location_exponent, next_available_complex_++);
	}


	void SLPCompiler::Visit(node::PowerOperator const& n){
		auto as_ptr = std::dynamic_pointer_cast<node::PowerOperator const>(n.shared_from_this());
		//get location of base and power then add instruction

		const auto& base = n.GetBase();
		const auto& exponent = n.GetExponent();

		if (this->locations_encountered_nodes_.find(base) == this->locations_encountered_nodes_.end())
			base->Accept(*this);

		if (this->locations_encountered_nodes_.find(exponent) == this->locations_encountered_nodes_.end())
			exponent->Accept(*this);

		auto loc_base = locations_encountered_nodes_[base];
		auto loc_exponent = locations_encountered_nodes_[exponent];

		this->locations_encountered_nodes_[as_ptr] =  next_available_complex_;
		program_under_construction_.AddInstruction(Power, loc_base, loc_exponent, next_available_complex_++);



	}

	void SLPCompiler::Visit(node::ExpOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::ExpOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Exp,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::LogOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::LogOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Log,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::NegateOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::NegateOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Negate,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::SqrtOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::SqrtOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Sqrt,location_operand, next_available_complex_++);
	}


	// the trig operators
	void SLPCompiler::Visit(node::SinOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::SinOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Sin,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::ArcSinOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::ArcSinOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Asin,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::CosOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::CosOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Cos,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::ArcCosOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::ArcCosOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Acos,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::TanOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::TanOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Tan,location_operand, next_available_complex_++);
	}

	void SLPCompiler::Visit(node::ArcTanOperator const& n){


		auto operand = n.Operand();
		if (this->locations_encountered_nodes_.find(operand) == this->locations_encountered_nodes_.end())
			operand->Accept(*this);

		auto location_operand = locations_encountered_nodes_[operand];
		this->locations_encountered_nodes_[std::dynamic_pointer_cast<node::ArcTanOperator const>(n.shared_from_this())] =  next_available_complex_;
		program_under_construction_.AddInstruction(Atan,location_operand, next_available_complex_++);
	}












	template <typename SourceT>
	SLP SLPCompiler::Compile(SourceT const& sys){
		this->Clear();

		// ThreadPrecision (thread-local), not DefaultPrecision (global): SLPs are
		// (re)compiled lazily during Eval, which may run on a std::thread worker
		// whose precision was set via SetThreadPrecision.  The global can be stale
		// (e.g. still at the boost default) on MPI worker ranks.
		const unsigned prec = ThreadPrecision();

		// deal with variables


			// 1. ADD VARIABLES

		program_under_construction_.input_locations_.Variables = next_available_complex_;

		auto variable_ordering = sys.VariableOrdering();
		for (auto v: variable_ordering){
			locations_encountered_nodes_[ v ] = next_available_complex_++;
		}
		program_under_construction_.number_of_.Variables = variable_ordering.size();

			// deal with path variable
		if (sys.HavePathVariable())
		{
				// do this action only if the system has a path variable defined
			program_under_construction_.input_locations_.Time = next_available_complex_;
			locations_encountered_nodes_[ sys.GetPathVariable() ] = next_available_complex_++;
			program_under_construction_.has_path_variable_ = true;
		}



			// 3. ADD FUNCTIONS AND DERIVATIVES (we omit the patches).
			//
			// Each output is a bare expression root: the natural functions' entry expressions,
			// then the space derivatives, then the time derivatives.  We reserve a contiguous
			// output slot for every output, then visit each root (computing its value into its
			// own slot) and emit an Assign copying that value into the reserved output slot.  The
			// compiler marks entry points itself via this explicit output list, rather than
			// relying on a Function wrapper node (ADR-0027).

		std::vector<std::shared_ptr<node::Node>> function_roots;
		for (auto const& f : sys.GetNaturalFunctions())
			function_roots.push_back(f);

		auto ds_dx = sys.GetSpaceDerivatives();
		std::vector<std::shared_ptr<node::Node>> ds_dt;
		if (sys.HavePathVariable())
			for (auto const& d : sys.GetTimeDerivatives())
				ds_dt.push_back(d);

		// reserve the output slots, contiguously, in the order [functions | jacobian | timederiv]
		program_under_construction_.number_of_.Functions = function_roots.size();
		program_under_construction_.output_locations_.Functions = next_available_complex_;
		std::vector<size_t> function_output_slots;
		for (size_t i = 0; i < function_roots.size(); ++i)
			function_output_slots.push_back(next_available_complex_++);

		program_under_construction_.number_of_.Jacobian = ds_dx.size();
		program_under_construction_.output_locations_.Jacobian = next_available_complex_;
		std::vector<size_t> jacobian_output_slots;
		for (size_t i = 0; i < ds_dx.size(); ++i)
			jacobian_output_slots.push_back(next_available_complex_++);

		std::vector<size_t> time_deriv_output_slots;
		if (sys.HavePathVariable()) {
			program_under_construction_.number_of_.TimeDeriv = ds_dt.size();
			program_under_construction_.output_locations_.TimeDeriv = next_available_complex_;
			for (size_t i = 0; i < ds_dt.size(); ++i)
				time_deriv_output_slots.push_back(next_available_complex_++);
		}

		// visit each output root and copy its value into the reserved output slot.  A shared root
		// is visited once (CSE), but every output position gets its own Assign.
		auto wire_outputs = [&](auto const& roots, std::vector<size_t> const& out_slots) {
			for (size_t i = 0; i < roots.size(); ++i) {
				auto const& r = roots[i];
				if (this->locations_encountered_nodes_.find(r) == this->locations_encountered_nodes_.end())
					r->Accept(*this);
				program_under_construction_.AddInstruction(Assign, this->locations_encountered_nodes_[r], out_slots[i]);
			}
		};

		wire_outputs(function_roots, function_output_slots);
		wire_outputs(ds_dx, jacobian_output_slots);
		if (sys.HavePathVariable())
			wire_outputs(ds_dt, time_deriv_output_slots);


		// the program's total slot count is the number of memory slots allocated during compilation
		program_under_construction_.num_slots_ = next_available_complex_;

		// NumType per slot (ADR-0034): infer which slots are real vs complex from the tree/tape.
		program_under_construction_.ComputeSlotNumTypes();

		// Split the tape into a frozen (constants-only) prologue and a live segment, so a
		// point-only re-evaluation can skip recomputing the constants (ADR-0027).  This is a pure
		// program operation (no memory needed).
		program_under_construction_.PartitionInstructions();


		// Wrap the now-immutable program in a facade and set up a per-thread memory for it.
		SLP result;
		result.program_ = std::make_shared<const SLPProgram>(std::move(program_under_construction_));
		result.memory_.precision_ = prec;
		result.SetupMemory();

		return result;
	}

	// Explicit instantiations: compile from a whole System (classic / slp_test path) and from a
	// PolynomialBlock (the fold -- the block compiles its own SLP in PolynomialBlock::Differentiate).
	template SLP SLPCompiler::Compile<System>(System const&);
	template SLP SLPCompiler::Compile<blocks::PolynomialBlock>(blocks::PolynomialBlock const&);

	void SLPCompiler::Clear(){
		next_available_complex_ = 0;
		next_available_int_ = 0;

		locations_encountered_nodes_.clear();
		program_under_construction_ = SLPProgram();
	}

}
