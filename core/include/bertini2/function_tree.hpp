//This file is part of Bertini 2.
//
//function_tree.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//function_tree.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with function_tree.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, University of Wisconsin Eau Claire
// jeb collins, west texas a&m

/**
\file function_tree.hpp 

\brief Collects the various header files which define the Bertini2 function tree node types.
*/


#ifndef BERTINI_FUNCTION_TREE_HPP
#define BERTINI_FUNCTION_TREE_HPP

#include "bertini2/function_tree/node.hpp"

#include "bertini2/function_tree/operators/operator.hpp"

#include "bertini2/function_tree/operators/arithmetic.hpp"
#include "bertini2/function_tree/operators/trig.hpp"
#include "bertini2/function_tree/symbols/linear_product.hpp"


#include "bertini2/function_tree/symbols/symbol.hpp"
#include "bertini2/function_tree/symbols/variable.hpp"
#include "bertini2/function_tree/symbols/number.hpp"
#include "bertini2/function_tree/symbols/special_number.hpp"

#include "bertini2/function_tree/roots/function.hpp"
#include "bertini2/function_tree/roots/jacobian.hpp"

#include "bertini2/function_tree/simplify.hpp"





BOOST_CLASS_EXPORT_KEY(bertini::node::Node)

BOOST_SERIALIZATION_ASSUME_ABSTRACT(bertini::node::Node)


BOOST_CLASS_EXPORT_KEY(bertini::node::Handle)

BOOST_CLASS_EXPORT_KEY(bertini::node::Symbol)
BOOST_CLASS_EXPORT_KEY(bertini::node::NamedSymbol)
BOOST_CLASS_EXPORT_KEY(bertini::node::Number)


BOOST_SERIALIZATION_ASSUME_ABSTRACT(bertini::node::Symbol)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(bertini::node::NamedSymbol)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(bertini::node::Number)


BOOST_CLASS_EXPORT_KEY(bertini::node::Variable)
BOOST_CLASS_EXPORT_KEY(bertini::node::Differential)

BOOST_CLASS_EXPORT_KEY(bertini::node::Float)
BOOST_CLASS_EXPORT_KEY(bertini::node::Integer)
BOOST_CLASS_EXPORT_KEY(bertini::node::Rational)

BOOST_CLASS_EXPORT_KEY(bertini::node::special_number::Pi)
BOOST_CLASS_EXPORT_KEY(bertini::node::special_number::E)

BOOST_CLASS_EXPORT_KEY(bertini::node::Function)
BOOST_CLASS_EXPORT_KEY(bertini::node::Jacobian)

BOOST_CLASS_EXPORT_KEY(bertini::node::SinOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::ArcSinOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::CosOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::ArcCosOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::TanOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::ArcTanOperator)


BOOST_CLASS_EXPORT_KEY(bertini::node::SumOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::NegateOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::MultOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::PowerOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::IntegerPowerOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::SqrtOperator)
BOOST_CLASS_EXPORT_KEY(bertini::node::ExpOperator)






BOOST_CLASS_TRACKING(bertini::node::Node, track_always)

BOOST_CLASS_TRACKING(bertini::node::Symbol, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::NamedSymbol, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::Number, boost::serialization::track_always)


BOOST_CLASS_TRACKING(bertini::node::Variable, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::Differential, boost::serialization::track_always)

BOOST_CLASS_TRACKING(bertini::node::Float, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::Integer, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::Rational, boost::serialization::track_always)

BOOST_CLASS_TRACKING(bertini::node::special_number::Pi, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::special_number::E, boost::serialization::track_always)

BOOST_CLASS_TRACKING(bertini::node::Function, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::Jacobian, boost::serialization::track_always)

BOOST_CLASS_TRACKING(bertini::node::SinOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::ArcSinOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::CosOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::ArcCosOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::TanOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::ArcTanOperator, boost::serialization::track_always)


BOOST_CLASS_TRACKING(bertini::node::SumOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::NegateOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::MultOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::PowerOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::IntegerPowerOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::SqrtOperator, boost::serialization::track_always)
BOOST_CLASS_TRACKING(bertini::node::ExpOperator, boost::serialization::track_always)





#endif




