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
// Copyright(C) 2016 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire


/**
\file include/bertini2/function_tree/forward_declares.hpp

\brief Forward declarations for types in the node:: namespace
*/

#pragma once

namespace bertini {
	namespace node{
		class Node;
	}
	
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
		class NegateOperator;
		class SqrtOperator;

		class LinearProduct;
		class DiffLinear;

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

	namespace node{
		namespace special_number{
			class Pi;
			class E;
		}
	}


	namespace node{
		template <typename Archive>
		void register_derived_node_types(Archive& ar)
		{
			ar.template register_type<bertini::node::Variable>();
			ar.template register_type<bertini::node::Integer>();
			ar.template register_type<bertini::node::Float>();
			ar.template register_type<bertini::node::Rational>();
			ar.template register_type<bertini::node::Function>();
			ar.template register_type<bertini::node::Jacobian>();
			ar.template register_type<bertini::node::Differential>();
			// ar.template register_type<bertini::node::Operator>(); // abstract type
			// ar.template register_type<bertini::node::UnaryOperator>(); // abstract type
			// ar.template register_type<bertini::node::NaryOperator>(); // abstract type
			ar.template register_type<bertini::node::SumOperator>();
			ar.template register_type<bertini::node::MultOperator>();
			ar.template register_type<bertini::node::IntegerPowerOperator>();
			ar.template register_type<bertini::node::PowerOperator>();
			ar.template register_type<bertini::node::ExpOperator>();
			ar.template register_type<bertini::node::LogOperator>();
			ar.template register_type<bertini::node::NegateOperator>();
			ar.template register_type<bertini::node::SqrtOperator>();
			ar.template register_type<bertini::node::LinearProduct>();
			// ar.template register_type<bertini::node::TrigOperator>(); // abstract type
			ar.template register_type<bertini::node::SinOperator>();
			ar.template register_type<bertini::node::ArcSinOperator>();
			ar.template register_type<bertini::node::CosOperator>();
			ar.template register_type<bertini::node::ArcCosOperator>();
			ar.template register_type<bertini::node::TanOperator>();
			ar.template register_type<bertini::node::ArcTanOperator>();
			ar.template register_type<bertini::node::special_number::Pi>();
			ar.template register_type<bertini::node::special_number::E>();
		}
	}
}// namespace bertini





