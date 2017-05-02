//This file is part of Bertini 2.
//
//src/function_tree/symbols/number.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/function_tree/symbols/number.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/function_tree/symbols/number.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of notre dame
// Jeb Collins, West Texas A&M


#include "function_tree/symbols/number.hpp"




namespace bertini{
	namespace node{
		using ::pow;
		
void Number::Reset() const
{
	// ResetStoredValues();
}


void Number::precision(unsigned int prec) const
{
	auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
	val_pair.first.precision(prec);
	val_pair.second = false; // false indicates to re-evaluate 
}


std::shared_ptr<Node> Number::Differentiate(std::shared_ptr<Variable> const& v) const
{
	return MakeInteger(0);
}

///////////////////
//
//  INTEGERS
//
////////////////////////


void Integer::print(std::ostream & target) const
{
	target << true_value_;
}

// Return value of constant
dbl Integer::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return dbl(double(true_value_),0);
}

void Integer::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = dbl(double(true_value_),0);
}


mpfr Integer::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return mpfr(true_value_,0);
}

void Integer::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = true_value_;
}


/////////////
//
//  Floats
//
////////////////


void Float::print(std::ostream & target) const
{
	target << highest_precision_value_;
}

// Return value of constant
dbl Float::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return dbl(highest_precision_value_);
}

void Float::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = dbl(highest_precision_value_);
}


mpfr Float::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return mpfr(highest_precision_value_);
}

void Float::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = mpfr(highest_precision_value_);
}


//
//  Rational
//



void Rational::print(std::ostream & target) const
{
	target << "(" << true_value_real_ << "," << true_value_imag_ << ")";
}

// Return value of constant
dbl Rational::FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const
{
	return dbl(double(true_value_real_),double(true_value_imag_));
}

void Rational::FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = dbl(double(true_value_real_),double(true_value_imag_));
}


mpfr Rational::FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const
{
	return mpfr(boost::multiprecision::mpfr_float(true_value_real_),boost::multiprecision::mpfr_float(true_value_imag_));
}

void Rational::FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const
{
	evaluation_value = mpfr(boost::multiprecision::mpfr_float(true_value_real_),boost::multiprecision::mpfr_float(true_value_imag_));
}


	} // re: namespace node
} // re: namespace bertini
