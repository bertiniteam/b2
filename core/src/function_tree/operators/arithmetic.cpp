//This file is part of Bertini 2.
//
//arithmetic.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//arithmetic.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with arithmetic.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire
// Jeb Collins, West Texas A&M


#include "bertini2/function_tree/operators/arithmetic.hpp"




namespace bertini{
	namespace node{
		using ::pow;
		
///////////////////////
//
//  Sum Operator definitions
//
//////////////////////





namespace{
	// print an operand, wrapping in parentheses only when its precedence is
	// too low for the position it occupies
	void PrintOperand(std::ostream& target, std::shared_ptr<Node> const& n, bool needs_parens)
	{
		if (needs_parens)
			target << "(";
		n->print(target);
		if (needs_parens)
			target << ")";
	}
}


std::shared_ptr<Node> SimplifiedNegate(std::shared_ptr<Node> const& n)
{
	if (n->IsLiteralZero())
		return n;
	return NegateOperator::Make(n);
}


std::shared_ptr<Node> SimplifiedSum(std::vector<std::pair<std::shared_ptr<Node>, bool>> const& terms)
{
	std::vector<std::pair<std::shared_ptr<Node>, bool>> remaining;
	remaining.reserve(terms.size());
	for (auto const& t : terms)
		if (!t.first->IsLiteralZero())
			remaining.push_back(t);

	if (remaining.empty())
		return Zero();
	if (remaining.size() == 1)
		return remaining[0].second ? remaining[0].first : SimplifiedNegate(remaining[0].first);

	return SumOperator::Make(remaining);  // build the complete sum, then intern once
}


namespace{
	// flatten nested Mult factors into one list so constants can meet and fold
	// (3*(2*x) -> {3,2,x} -> 6*x).  reads the nested node's operands; never
	// modifies it.  a divided nested Mult inverts its flags: x/(a/b) = x/a*b.
	void FlattenFactor(std::vector<std::pair<std::shared_ptr<Node>, bool>>& out,
	                   std::shared_ptr<Node> const& n, bool mult)
	{
		if (auto as_mult = std::dynamic_pointer_cast<MultOperator>(n))
		{
			auto const& ops = as_mult->Operands();
			auto const& flags = as_mult->GetMultOrDiv();
			for (size_t ii = 0; ii < ops.size(); ++ii)
				FlattenFactor(out, ops[ii], mult ? flags[ii] : !flags[ii]);
			return;
		}
		out.emplace_back(n, mult);
	}
}

std::shared_ptr<Node> SimplifiedMult(std::vector<std::pair<std::shared_ptr<Node>, bool>> const& factors_in)
{
	std::vector<std::pair<std::shared_ptr<Node>, bool>> factors;
	factors.reserve(factors_in.size());
	for (auto const& f : factors_in)
		FlattenFactor(factors, f.first, f.second);

	mpq_rational constant(1);
	bool have_constant = false;
	std::vector<std::pair<std::shared_ptr<Node>, bool>> remaining;
	remaining.reserve(factors.size());

	for (auto const& f : factors)
	{
		auto const& n = f.first;
		const bool mult = f.second;

		if (n->IsLiteralZero())
		{
			if (mult)
				return Zero();
			remaining.push_back(f);  // division by a literal zero stays visible
			continue;
		}
		if (n->IsLiteralOne())
			continue;

		// fold exact constants together; Floats are deliberately NOT folded
		if (auto as_int = std::dynamic_pointer_cast<Integer>(n))
		{
			mpq_rational val(as_int->GetValue());
			if (mult)
				constant *= val;
			else
				constant /= val;
			have_constant = true;
			continue;
		}
		if (auto as_rat = std::dynamic_pointer_cast<Rational>(n))
		{
			if (as_rat->GetValueImag() == 0)
			{
				if (mult)
					constant *= as_rat->GetValueReal();
				else
					constant /= as_rat->GetValueReal();
				have_constant = true;
				continue;
			}
		}
		remaining.push_back(f);
	}

	std::shared_ptr<Node> constant_node = nullptr;
	if (have_constant && constant != 1)
	{
		if (denominator(constant) == 1)
			constant_node = Integer::Make(numerator(constant));
		else
			constant_node = Rational::Make(constant, mpq_rational(0));
	}

	if (remaining.empty())
		return constant_node ? constant_node : std::shared_ptr<Node>(One());

	std::vector<std::pair<std::shared_ptr<Node>, bool>> finals;
	if (constant_node)
		finals.emplace_back(constant_node, true);  // canonical order: constant first
	finals.insert(finals.end(), remaining.begin(), remaining.end());

	if (finals.size() == 1 && finals[0].second)
		return finals[0].first;

	// for a leading divisor, materialize the canonical '1/...' form
	if (!finals[0].second)
		finals.insert(finals.begin(), {One(), true});

	return MultOperator::Make(finals);  // build the complete product, then intern once
}

// ---- functional Simplified() (non-mutating successors to Eliminate*/ReduceDepth) ----
// Each recurses on children (which return fresh simplified subtrees, sharing what they
// didn't change) and re-assembles through the Simplified* factories, so literal zeros/ones
// vanish, exact constants fold, and nested same-type operators flatten.

std::shared_ptr<Node> SumOperator::Simplified() const
{
	std::vector<std::pair<std::shared_ptr<Node>, bool>> terms;
	terms.reserve(operands_.size());
	for (size_t ii = 0; ii < operands_.size(); ++ii)
		terms.emplace_back(operands_[ii]->Simplified(), signs_[ii]);
	return SimplifiedSum(terms);
}

std::shared_ptr<Node> MultOperator::Simplified() const
{
	std::vector<std::pair<std::shared_ptr<Node>, bool>> factors;
	factors.reserve(operands_.size());
	for (size_t ii = 0; ii < operands_.size(); ++ii)
		factors.emplace_back(operands_[ii]->Simplified(), mult_or_div_[ii]);
	return SimplifiedMult(factors);
}

std::shared_ptr<Node> NegateOperator::Simplified() const
{
	return SimplifiedNegate(operand_->Simplified());
}

std::shared_ptr<Node> PowerOperator::Simplified() const
{
	auto base_s = base_->Simplified();
	auto exp_s  = exponent_->Simplified();
	if (exp_s->IsLiteralZero()) return Integer::Make(1);   // x^0 -> 1
	if (exp_s->IsLiteralOne())  return base_s;             // x^1 -> x
	return PowerOperator::Make(base_s, exp_s);
}

std::shared_ptr<Node> IntegerPowerOperator::Simplified() const
{
	if (exponent_ == 0) return Integer::Make(1);           // x^0 -> 1
	auto op_s = operand_->Simplified();
	if (exponent_ == 1) return op_s;                       // x^1 -> x
	return IntegerPowerOperator::Make(op_s, exponent_);
}

std::shared_ptr<Node> SqrtOperator::Simplified() const
{
	return SqrtOperator::Make(operand_->Simplified());
}

std::shared_ptr<Node> ExpOperator::Simplified() const
{
	return ExpOperator::Make(operand_->Simplified());
}

std::shared_ptr<Node> LogOperator::Simplified() const
{
	return LogOperator::Make(operand_->Simplified());
}

void SumOperator::print(std::ostream & target) const
{
	// Print a positive term first.  Canonical ordering (Rung 3c) sorts by degree, which can put
	// a subtracted term ahead of a constant ("1-t" -> operands t,1), and leading with that minus
	// reads as an extra negation ("-t+1").  Leading with a '+' term gives the natural, fewer-ops
	// form ("1-t", "x-(y+z)").  (Print order only; the canonical operand order is unchanged.  An
	// all-negative sum still leads with '-', e.g. "-x-y".)
	size_t lead = 0;
	for (size_t ii = 0; ii < operands_.size(); ++ii)
		if (signs_[ii]) { lead = ii; break; }

	auto print_one = [&](size_t ii, bool is_lead)
	{
		const bool plus = signs_[ii];
		if (is_lead)
		{
			if (!plus)
				target << "-";
		}
		else
			target << (plus ? "+" : "-");

		const auto prec = operands_[ii]->Precedence();
		// after '-', wrap sums (grouping) and anything printing a leading '-'
		// (avoids "--"); after '+' or in the lead, wrap only leading-'-' printers
		bool needs_parens;
		if (is_lead)
			needs_parens = plus ? false : (prec <= PrecNegate);
		else
			needs_parens = plus ? (prec == PrecNegate) : (prec <= PrecNegate);
		PrintOperand(target, operands_[ii], needs_parens);
	};

	print_one(lead, true);
	for (size_t ii = 0; ii < operands_.size(); ++ii)
		if (ii != lead)
			print_one(ii, false);
}


std::shared_ptr<Node> SumOperator::Differentiate(std::shared_ptr<Variable> const& v) const
{
	std::vector<std::pair<std::shared_ptr<Node>, bool>> terms;
	terms.reserve(operands_.size());
	for (size_t ii = 0; ii < operands_.size(); ++ii)
	{
		if (std::dynamic_pointer_cast<Number>(operands_[ii]))
			continue;  // constants differentiate to 0; don't even build it

		terms.emplace_back(operands_[ii]->Differentiate(v), signs_[ii]);
	}
	return SimplifiedSum(terms);
}

int SumOperator::Degree(std::shared_ptr<Variable> const& v) const
{
	int deg = 0;
	
	for (auto iter: operands_)
	{
		auto curr_deg = iter->Degree(v);
		if (curr_deg<0)
			return curr_deg;
		
		deg = std::max(deg, curr_deg);
	}
	return deg;
}

int SumOperator::Degree(VariableGroup const& vars) const 
{
	auto deg = 0;
	for (auto iter = operands_.begin(); iter!=operands_.end(); iter++)
	{
		auto term_degree = (*iter)->Degree(vars);
		if (term_degree<0)
			return term_degree;

		deg = std::max(deg, term_degree);
		
	}
	
	return deg;
}


std::vector<int> SumOperator::MultiDegree(VariableGroup const& vars) const
{
	std::vector<int> deg(vars.size(),0);
	for (auto iter : operands_)
	{
		auto term_deg = iter->MultiDegree(vars);

		for (auto iter = term_deg.begin(); iter!= term_deg.end(); ++iter)
		{
			*(deg.begin()+(iter-term_deg.begin())) = std::max(*(deg.begin()+(iter-term_deg.begin())), *iter);
		}
	}
	return deg;
}

std::shared_ptr<Node> SumOperator::Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const
{
	// homogenize each summand functionally (fresh subtrees), measure degrees, then pad
	// degree-deficient summands with powers of homvar -- all into a freshly-built sum.
	// the input is never touched, and the throw below happens before anything is built,
	// so a non-polynomial term can't leave a half-homogenized tree behind.
	std::vector<std::shared_ptr<Node>> new_ops;
	new_ops.reserve(operands_.size());
	for (auto const& op : operands_)
		new_ops.push_back(op->Homogenized(vars, homvar));

	int maxdegree = 0;
	std::vector<int> term_degrees;
	term_degrees.reserve(new_ops.size());
	for (auto const& op : new_ops)
	{
		auto local_degree = op->Degree(vars);
		if (local_degree < 0)
			throw std::runtime_error("asking for homogenization on non-polynomial node");
		term_degrees.push_back(local_degree);
		maxdegree = std::max(maxdegree, local_degree);
	}

	for (size_t ii = 0; ii < new_ops.size(); ++ii)
	{
		auto degree_deficiency = maxdegree - term_degrees[ii];
		if (degree_deficiency == 1)
			new_ops[ii] = MultOperator::Make(homvar, new_ops[ii]);
		else if (degree_deficiency > 1)
			new_ops[ii] = MultOperator::Make(
				IntegerPowerOperator::Make(std::static_pointer_cast<Node>(homvar), degree_deficiency),
				new_ops[ii]);
	}

	std::vector<std::pair<std::shared_ptr<Node>, bool>> terms;
	terms.reserve(new_ops.size());
	for (size_t ii = 0; ii < new_ops.size(); ++ii)
		terms.emplace_back(new_ops[ii], signs_[ii]);
	return SumOperator::Make(terms);  // complete sum, interned once
}


bool SumOperator::IsHomogeneous(std::shared_ptr<Variable> const& v) const
{
	
	for (auto iter : operands_)
	{
		if (!iter->IsHomogeneous(v))
			return false;
	}

	// the only hope this has of being homogeneous, is that each factor is homogeneous
	int deg;

	deg = (*(operands_.begin()))->Degree(v) ;

	if (deg < 0) 
		return false;

	for (auto iter = operands_.begin()+1; iter!= operands_.end(); iter++)
	{
		auto local_degree = (*iter)->Degree(v);
		if (local_degree!=deg)
			return false;
	}

	return true;
}

bool SumOperator::IsHomogeneous(VariableGroup const& v) const
{
	
	for (auto iter : operands_)
	{
		if (!iter->IsHomogeneous(v))
			return false;
	}

	// the only hope this has of being homogeneous, is that each factor is homogeneous
	int deg;

	deg = (*(operands_.begin()))->Degree(v) ;

	if (deg < 0) 
		return false;

	for (auto iter = operands_.begin()+1; iter!= operands_.end(); iter++)
	{
		auto local_degree = (*iter)->Degree(v);
		if (local_degree!=deg)
			return false;
	}

	return true;
}




dbl SumOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	dbl retval;
	this->FreshEval_d(retval, diff_variable);
	return retval;
}
	
void SumOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = dbl(0);
	for(size_t ii = 0; ii < operands_.size(); ++ii)
	{
		if(signs_[ii])
		{
			operands_[ii]->EvalInPlace<dbl>(temp_d_, diff_variable);
			evaluation_value += temp_d_;
		}
		else
		{
			operands_[ii]->EvalInPlace<dbl>(temp_d_, diff_variable);
			evaluation_value -= temp_d_;
		}
	}
}

	
	

mpfr_complex SumOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	mpfr_complex retval;
	this->FreshEval_mp(retval, diff_variable);
	return retval;
}

void SumOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	if (signs_[0])
		operands_[0]->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
	else
	{
		operands_[0]->EvalInPlace<mpfr_complex>(temp_mp_, diff_variable);
		evaluation_value = -temp_mp_;
	}

	for(size_t ii = 1; ii < operands_.size(); ++ii)
	{
		if(signs_[ii])
		{
			operands_[ii]->EvalInPlace<mpfr_complex>(temp_mp_, diff_variable);
			evaluation_value += temp_mp_;
		}
		else
		{
			operands_[ii]->EvalInPlace<mpfr_complex>(temp_mp_, diff_variable);
			evaluation_value -= temp_mp_;
		}
	}
	
}


















//////////////////////
//
//  Negate operator definitions
//
////////////////////////

void NegateOperator::print(std::ostream & target) const
{
	target << "-";
	PrintOperand(target, operand_, operand_->Precedence() <= PrecNegate);
}

std::shared_ptr<Node> NegateOperator::Differentiate(std::shared_ptr<Variable> const& v) const
{
	return SimplifiedNegate(operand_->Differentiate(v));
}

dbl NegateOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return -(operand_->Eval<dbl>(diff_variable));
}

void NegateOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
	evaluation_value = -evaluation_value;
}


mpfr_complex NegateOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return -operand_->Eval<mpfr_complex>(diff_variable);
}

void NegateOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
	evaluation_value = -evaluation_value;
}























///////////////////////
//
//  Mult Operator definitions
//
//////////////////////









void MultOperator::print(std::ostream & target) const
{
	for (size_t ii = 0; ii < operands_.size(); ++ii)
	{
		const bool mult = mult_or_div_[ii];
		if (ii == 0)
		{
			if (!mult)
				target << "1/";
		}
		else
			target << (mult ? "*" : "/");

		const auto prec = operands_[ii]->Precedence();
		// multiplied positions: wrap below-mult precedence (sums, leading-'-'
		// printers); divided positions: also wrap other mults (grouping)
		const bool needs_parens = mult ? (prec < PrecMult) : (prec <= PrecMult);
		PrintOperand(target, operands_[ii], needs_parens);
	}
}




std::shared_ptr<Node> MultOperator::Differentiate(std::shared_ptr<Variable> const& v) const
{
	std::vector<std::pair<std::shared_ptr<Node>, bool>> sum_terms;
	// this loop implements the generic product rule, perhaps inefficiently.
	for (size_t ii = 0; ii < operands_.size(); ++ii)
	{
		auto local_derivative = operands_[ii]->Differentiate(v);
		if (local_derivative->IsLiteralZero())
			continue;

		// the product of the derivative with the remaining factors
		std::vector<std::pair<std::shared_ptr<Node>, bool>> factors;
		factors.reserve(operands_.size() + 1);
		factors.emplace_back(local_derivative, true);
		for (size_t jj = 0; jj < operands_.size(); ++jj)
			if (jj != ii)
				factors.emplace_back(operands_[jj], mult_or_div_[jj]);

		// if is division, need this for the quotient rule
		if ( !(mult_or_div_[ii]) )
			factors.emplace_back(pow(operands_[ii],2), false); // draw a line and square below

		// a divided factor's term enters subtracted (quotient rule)
		sum_terms.emplace_back(SimplifiedMult(factors), mult_or_div_[ii]);
	} // re: for ii

	return SimplifiedSum(sum_terms);
}

int MultOperator::Degree(std::shared_ptr<Variable> const& v) const
{
	int deg = 0;
	for (auto iter = operands_.begin(); iter!= operands_.end(); iter++)
	{
		// if the operand node is a differential coming from another variable, then this degree should be 0, end of story.
		
		auto factor_deg = (*iter)->Degree(v);
		
		auto is_it_a_differential = std::dynamic_pointer_cast<Differential>(*iter);
		if (is_it_a_differential)
			if (is_it_a_differential->GetVariable()!=v)
			{
				return 0;
			}
		
		
		
		
		
		if (factor_deg<0)
			return factor_deg;
		else if (factor_deg!=0 && !*(mult_or_div_.begin() + (iter-operands_.begin()) ) )
			return -1;
		else
			deg+=factor_deg;
	}
	return deg;
}


int MultOperator::Degree(VariableGroup const& vars) const 
{
	
	auto deg = 0;

	for (auto iter = operands_.begin(); iter!=operands_.end(); iter++)
	{
		auto factor_deg = (*iter)->Degree(vars);  

		if (factor_deg<0)
			return factor_deg;
		else if (factor_deg!=0 && !*(mult_or_div_.begin() + (iter-operands_.begin()) ) )
			return -1;
		else
			deg+=factor_deg;
	}

	
	return deg;
}

std::vector<int> MultOperator::MultiDegree(VariableGroup const& vars) const
{
	std::vector<int> deg(vars.size(),0);
	for (auto iter : operands_)
	{
		auto term_deg = iter->MultiDegree(vars);

		for (auto iter = term_deg.begin(); iter!= term_deg.end(); ++iter)
		{
			*(deg.begin()+(iter-term_deg.begin())) += *iter;
		}
	}
	return deg;
}



std::shared_ptr<Node> MultOperator::Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const
{
	// product of homogenized factors, preserving the multiply/divide flags (no simplification).
	std::vector<std::shared_ptr<Node>> ops;
	ops.reserve(operands_.size());
	for (auto const& op : operands_)
		ops.push_back(op->Homogenized(vars, homvar));

	std::vector<std::pair<std::shared_ptr<Node>, bool>> factors;
	factors.reserve(ops.size() + 1);
	if (!mult_or_div_[0])
		factors.emplace_back(One(), true);   // leading divisor -> canonical 1/...
	for (size_t ii = 0; ii < ops.size(); ++ii)
		factors.emplace_back(ops[ii], mult_or_div_[ii]);
	return MultOperator::Make(factors);  // complete product, interned once
}


bool MultOperator::IsHomogeneous(std::shared_ptr<Variable> const& v) const
{
	// the only hope this has of being homogeneous, is that each factor is homogeneous
	for (auto iter : operands_)
	{
		if (! iter->IsHomogeneous(v))
		{
			return false;
		}
	}
	return true;
}

bool MultOperator::IsHomogeneous(VariableGroup const& v) const
{
	// the only hope this has of being homogeneous, is that each factor is homogeneous
	for (auto iter : operands_)
	{
		if (! iter->IsHomogeneous(v))
		{
			return false;
		}
	}
	return true;
}

dbl MultOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	dbl retval;
	this->FreshEval_d(retval, diff_variable);
	return retval;
}

void MultOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = dbl(1);
	for(size_t ii = 0; ii < operands_.size(); ++ii)
	{
		if(mult_or_div_[ii])
		{
			operands_[ii]->EvalInPlace<dbl>(temp_d_, diff_variable);
			evaluation_value *= temp_d_;
		}
		else
		{
			operands_[ii]->EvalInPlace<dbl>(temp_d_, diff_variable);
			evaluation_value /= temp_d_;
		}
	}
	
}


mpfr_complex MultOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	mpfr_complex retval;
	this->FreshEval_mp(retval, diff_variable);
	return retval;
}

void MultOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	if (mult_or_div_[0])
		operands_[0]->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
	else
	{
		operands_[0]->EvalInPlace<mpfr_complex>(temp_mp_, diff_variable);
		evaluation_value = static_cast<mpfr_float>(1)/temp_mp_;
	}

	for(size_t ii = 1; ii < operands_.size(); ++ii)
	{
		operands_[ii]->EvalInPlace<mpfr_complex>(temp_mp_, diff_variable);
		if(mult_or_div_[ii])
			evaluation_value *= temp_mp_;
		else
			evaluation_value /= temp_mp_;
	}
	
}



////////////
//
//  Power Operator definitions
//
/////////////////



void PowerOperator::Reset() const
{
	Node::ResetStoredValues();
	base_->Reset();
	exponent_->Reset();
}

void PowerOperator::print(std::ostream & target) const
{
	// '^' is right-associative and binds tightest, so wrap anything that is
	// not an atom -- including other powers, to keep x^y^z unambiguous
	PrintOperand(target, base_, base_->Precedence() <= PrecPower);
	target << "^";
	PrintOperand(target, exponent_, exponent_->Precedence() <= PrecPower);
}


std::shared_ptr<Node> PowerOperator::Differentiate(std::shared_ptr<Variable> const& v) const
{
	auto exp_minus_one = exponent_-1;
	return SimplifiedMult({
		{base_->Differentiate(v), true},
		{exponent_, true},
		{PowerOperator::Make(base_, exp_minus_one), true}
	});
}


int PowerOperator::Degree(std::shared_ptr<Variable> const& v) const
{
	
	auto base_deg = base_->Degree(v);
	auto exp_deg = exponent_->Degree(v);
	
	if (exp_deg==0)
	{
		auto exp_val = exponent_->Eval<dbl>();
		bool exp_is_int = false;
		
		if (fabs(imag(exp_val))< 10*std::numeric_limits<double>::epsilon()) // so a real thresholding step
			if (fabs(real(exp_val) - std::round(real(exp_val))) < 10*std::numeric_limits<double>::epsilon()) // then check that the real part is close to an integer
				exp_is_int = true;
		
		if (exp_is_int)
		{
			if (abs(exp_val-dbl(0.0))< 10*std::numeric_limits<double>::epsilon())
				return 0;
			else if (real(exp_val)<0)
				return -1;
			else  // positive integer.
			{
				if (base_deg<0)
					return -1;
				else
					return base_deg*static_cast<int>(std::round(real(exp_val)));
			}
			
		}
		else
		{
			if (base_deg==0)
				return 0;
			else
				return -1;
		}
		
	}
	else
	{
		// there may be an edge case here where the base is the number 0 or 1.  but that would be stupid, wouldn't it.
		return -1;
	}
}

int PowerOperator::Degree(VariableGroup const& vars) const
{
	auto multideg = MultiDegree(vars);
	auto deg = 0;
	std::for_each(multideg.begin(),multideg.end(),[&](int n){
					if (n < 0)
						deg = -1;
					else
						deg += n;
					});
	return deg;
}


std::vector<int> PowerOperator::MultiDegree(VariableGroup const& vars) const
{
	std::vector<int> deg(vars.size(),0);
	for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
	{
		*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
	}
	return deg;
}



std::shared_ptr<Node> PowerOperator::Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const
{
	if (exponent_->Degree(vars) == 0)
		return PowerOperator::Make(base_->Homogenized(vars, homvar), exponent_);
	// non-constant exponent -> non-polynomial; throw before building anything.
	throw std::runtime_error("asking for homogenization on non-polynomial node");
}

// ---- unary operators: rebuild the same op type from the homogenized operand ----
std::shared_ptr<Node> NegateOperator::Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const
{
	return NegateOperator::Make(operand_->Homogenized(vars, homvar));
}

std::shared_ptr<Node> IntegerPowerOperator::Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const
{
	return IntegerPowerOperator::Make(operand_->Homogenized(vars, homvar), exponent_);
}

std::shared_ptr<Node> SqrtOperator::Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const
{
	return SqrtOperator::Make(operand_->Homogenized(vars, homvar));
}

std::shared_ptr<Node> ExpOperator::Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const
{
	return ExpOperator::Make(operand_->Homogenized(vars, homvar));
}

std::shared_ptr<Node> LogOperator::Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const
{
	return LogOperator::Make(operand_->Homogenized(vars, homvar));
}

// ---- structural hash / equality (Rung 2) ----
// Order-sensitive; operands folded in by Hash() for the hash, compared by pointer for IsSame.

std::size_t SumOperator::HashImpl() const
{
	std::size_t h = typeid(SumOperator).hash_code();
	for (size_t ii = 0; ii < operands_.size(); ++ii)
	{
		Node::HashCombine(h, operands_[ii]->Hash());
		Node::HashCombine(h, signs_[ii] ? 1u : 0u);
	}
	return h;
}
bool SumOperator::IsSame(Node const& other) const
{
	auto o = dynamic_cast<SumOperator const*>(&other);
	if (!o || operands_.size() != o->operands_.size() || signs_ != o->signs_)
		return false;
	for (size_t ii = 0; ii < operands_.size(); ++ii)
		if (operands_[ii].get() != o->operands_[ii].get())
			return false;
	return true;
}

std::size_t MultOperator::HashImpl() const
{
	std::size_t h = typeid(MultOperator).hash_code();
	for (size_t ii = 0; ii < operands_.size(); ++ii)
	{
		Node::HashCombine(h, operands_[ii]->Hash());
		Node::HashCombine(h, mult_or_div_[ii] ? 1u : 0u);
	}
	return h;
}
bool MultOperator::IsSame(Node const& other) const
{
	auto o = dynamic_cast<MultOperator const*>(&other);
	if (!o || operands_.size() != o->operands_.size() || mult_or_div_ != o->mult_or_div_)
		return false;
	for (size_t ii = 0; ii < operands_.size(); ++ii)
		if (operands_[ii].get() != o->operands_[ii].get())
			return false;
	return true;
}

std::size_t PowerOperator::HashImpl() const
{
	std::size_t h = typeid(PowerOperator).hash_code();
	Node::HashCombine(h, base_->Hash());
	Node::HashCombine(h, exponent_->Hash());
	return h;
}
bool PowerOperator::IsSame(Node const& other) const
{
	auto o = dynamic_cast<PowerOperator const*>(&other);
	return o && base_.get() == o->base_.get() && exponent_.get() == o->exponent_.get();
}

std::size_t IntegerPowerOperator::HashImpl() const
{
	std::size_t h = typeid(IntegerPowerOperator).hash_code();
	Node::HashCombine(h, operand_->Hash());
	Node::HashCombine(h, std::hash<int>{}(exponent_));
	return h;
}
bool IntegerPowerOperator::IsSame(Node const& other) const
{
	auto o = dynamic_cast<IntegerPowerOperator const*>(&other);
	return o && exponent_ == o->exponent_ && operand_.get() == o->operand_.get();
}

bool PowerOperator::IsHomogeneous(std::shared_ptr<Variable> const& v) const
{
	// the only hope this has of being homogeneous, is that the degree of the exponent is 0 (it's constant), and that it's an integer
	if (exponent_->Degree(v)==0)
	{
		auto exp_val = exponent_->Eval<dbl>();
		if (fabs(imag(exp_val)) < 10*std::numeric_limits<double>::epsilon())
			if (fabs(std::round(real(exp_val)) - real(exp_val)) < 10*std::numeric_limits<double>::epsilon())
				if (real(exp_val) >=0 )
					return base_->IsHomogeneous(v);
	}
	return false;
}


bool PowerOperator::IsHomogeneous(VariableGroup const& v) const
{
	// the only hope this has of being homogeneous, is that the degree of the exponent is 0 (it's constant), and that it's an integer
	if (exponent_->Degree(v)==0)
	{
		auto exp_val = exponent_->Eval<dbl>();
		if (fabs(imag(exp_val)) < 10*std::numeric_limits<double>::epsilon())
			if (fabs(std::round(real(exp_val)) - real(exp_val)) < 10*std::numeric_limits<double>::epsilon())
				if (real(exp_val) >=0 )
					return base_->IsHomogeneous(v);
	}
	return false;
}

dbl PowerOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return std::pow( base_->Eval<dbl>(diff_variable), exponent_->Eval<dbl>());
}

void PowerOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	dbl temp_d;
	exponent_->EvalInPlace<dbl>(temp_d);
	base_->EvalInPlace<dbl>(evaluation_value, diff_variable);
	
	evaluation_value = std::pow(evaluation_value, temp_d);
}


mpfr_complex PowerOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return pow( base_->Eval<mpfr_complex>(diff_variable), exponent_->Eval<mpfr_complex>());
}

void PowerOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	mpfr_complex temp_mp;
	exponent_->EvalInPlace<mpfr_complex>(temp_mp);
	base_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
	
	evaluation_value = pow(evaluation_value, temp_mp);
}
















/////////////////
//
//  IntegerPowerOperator definitions
//
////////////////////

void IntegerPowerOperator::print(std::ostream & target) const
{
	PrintOperand(target, operand_, operand_->Precedence() <= PrecPower);
	if (exponent() < 0)
		target << "^(" << exponent() << ")";
	else
		target << "^" << exponent();
}



std::shared_ptr<Node> IntegerPowerOperator::Differentiate(std::shared_ptr<Variable> const& v) const
{
	if (exponent_==0)
		return Integer::Make(0);
	else if (exponent_==1)
		return operand_->Differentiate(v);
	else{
		std::shared_ptr<Node> power_part = (exponent_==2) ? operand_ : std::shared_ptr<Node>(IntegerPowerOperator::Make(operand_, exponent_-1));
		return SimplifiedMult({
			{Integer::Make(exponent_), true},
			{power_part, true},
			{operand_->Differentiate(v), true}
		});
	}
}


int IntegerPowerOperator::Degree(std::shared_ptr<Variable> const& v) const
{
	auto base_deg = operand_->Degree(v);
	if (base_deg<0)
		return base_deg;
	else
		return exponent_*base_deg;
	
}












//////////////
//
//  Square Root Operator definitions
//
/////////////////

void SqrtOperator::print(std::ostream & target) const
{
	target << "sqrt(";
	operand_->print(target);
	target << ")";
}



std::shared_ptr<Node> SqrtOperator::Differentiate(std::shared_ptr<Variable> const& v) const
{
	return SimplifiedMult({
		{Rational::Make(mpq_rational(1,2),0), true},
		{PowerOperator::Make(operand_, Rational::Make(mpq_rational(-1,2),0)), true},
		{operand_->Differentiate(v), true}
	});
}

int SqrtOperator::Degree(std::shared_ptr<Variable> const& v) const
{
	if (operand_->Degree(v)==0)
	{
		return 0;
	}
	else
	{
		return -1;
	}
}


dbl SqrtOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return sqrt(operand_->Eval<dbl>(diff_variable));
}

void SqrtOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
	evaluation_value = sqrt(evaluation_value);
}


mpfr_complex SqrtOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return sqrt(operand_->Eval<mpfr_complex>(diff_variable));
}

void SqrtOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
	evaluation_value = sqrt(evaluation_value);
}















///////////////
//
//  ExpOperator definitions
//
//////////////

void ExpOperator::print(std::ostream & target) const
{
	target << "exp(";
	operand_->print(target);
	target << ")";
}




std::shared_ptr<Node> ExpOperator::Differentiate(std::shared_ptr<Variable> const& v) const
{
	return SimplifiedMult({
		{exp(operand_), true},
		{operand_->Differentiate(v), true}
	});
}


int ExpOperator::Degree(std::shared_ptr<Variable> const& v) const
{
	if (operand_->Degree(v)==0)
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

dbl ExpOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return exp(operand_->Eval<dbl>(diff_variable));
}

void ExpOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
	evaluation_value = exp(evaluation_value);
}


mpfr_complex ExpOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return exp(operand_->Eval<mpfr_complex>(diff_variable));
}

void ExpOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
	evaluation_value = exp(evaluation_value);
}











///////////////
//
//  LogOperator definitions
//
//////////////

void LogOperator::print(std::ostream & target) const
{
	target << "log(";
	operand_->print(target);
	target << ")";
}



std::shared_ptr<Node> LogOperator::Differentiate(std::shared_ptr<Variable> const& v) const
{
	return SimplifiedMult({
		{operand_->Differentiate(v), true},
		{operand_, false}
	});
}


int LogOperator::Degree(std::shared_ptr<Variable> const& v) const
{
	if (operand_->Degree(v)==0)
	{
		return 0;
	}
	else
	{
		return -1;
	}
}

dbl LogOperator::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return log(operand_->Eval<dbl>(diff_variable));
}

void LogOperator::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	operand_->EvalInPlace<dbl>(evaluation_value, diff_variable);
	evaluation_value = log(evaluation_value);
}


mpfr_complex LogOperator::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return log(operand_->Eval<mpfr_complex>(diff_variable));
}

void LogOperator::FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	operand_->EvalInPlace<mpfr_complex>(evaluation_value, diff_variable);
	evaluation_value = log(evaluation_value);
}


	} // re: namespace node
} // re: namespace bertini
