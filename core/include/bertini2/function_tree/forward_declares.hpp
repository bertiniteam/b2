//This file is part of Bertini 2.
//
//include/bertini2/function_tree/forward_declares.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//include/bertini2/function_tree/forward_declares.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with include/bertini2/function_tree/forward_declares.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// Dani Brake
// University of Notre Dame


/**
\file include/bertini2/function_tree/forward_declares.hpp

\brief Forward declarations for types in the node:: namespace
*/

#pragma once

namespace bertini {

	namespace node{ 
		class Variable;
		class Integer;
		class Float;
		class Rational;
		class Function;
		class Jacobian;
		class Differential;
	}


	namespace node{ 
		class Operator;

		class UnaryOperator;
		class NaryOperator;
		
		class SumOperator;
		class MultOperator;
		class IntegerPowerOperator;
		class PowerOperator;
		class ExpOperator;
		class LogOperator;
	}

	namespace node{
		class TrigOperator;

		class SinOperator;
		class ArcSinOperator;
		class CosOperator;
		class ArcCosOperator;
		class TanOperator;
		class ArcTanOperator;
	}
}// namespace bertini
