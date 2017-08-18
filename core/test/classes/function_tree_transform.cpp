//This file is part of Bertini 2.
//
//b2/core/test/classes/function_tree_transform.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//b2/core/test/classes/function_tree_transform.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with b2/core/test/classes/function_tree_transform.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  Dani Brake
//  University of Wisconsin Eau Claire
//  ACMS
//  Fall 2017

#include <iostream>

#include <cstdlib>
#include <cmath>

#include "bertini2/function_tree.hpp"


#include <boost/spirit/include/qi.hpp>
#include <boost/test/unit_test.hpp>

#include "externs.hpp"

using Nd = std::shared_ptr<bertini::node::Node>;
using Nd = std::shared_ptr<bertini::node::Node>;

using bertini::MakeVariable;
using bertini::MakeInteger;

BOOST_AUTO_TEST_SUITE(function_tree)

BOOST_AUTO_TEST_SUITE(transform)

BOOST_AUTO_TEST_SUITE(eliminate_zeros)


BOOST_AUTO_TEST_CASE(cant_eliminate_self)
{
	auto zero = Nd(MakeInteger(0));

	auto num_eliminated = zero->EliminateZeros();

	BOOST_CHECK_EQUAL(num_eliminated, 0);
}

BOOST_AUTO_TEST_CASE(cant_eliminate_self2)
{
	auto zero = Nd(MakeInteger(0));
	auto n = zero;

	auto num_eliminated = n->EliminateZeros();

	BOOST_CHECK_EQUAL(num_eliminated, 0);
}


BOOST_AUTO_TEST_CASE(eliminate_single_zero_sum)
{
	auto zero = Nd(MakeInteger(0));

	auto n = zero+zero;

	auto num_eliminated = n->EliminateZeros();

	BOOST_CHECK_EQUAL(num_eliminated, 1);
}


BOOST_AUTO_TEST_CASE(eliminate_single_zero_sum2)
{
	auto zero = Nd(MakeInteger(0));

	auto n = zero+0;

	auto num_eliminated = n->EliminateZeros();

	BOOST_CHECK_EQUAL(num_eliminated, 1);
}


BOOST_AUTO_TEST_CASE(eliminate_single_zero_sum3)
{
	auto zero = Nd(MakeInteger(0));

	auto n = 0+zero;

	auto num_eliminated = n->EliminateZeros();

	BOOST_CHECK_EQUAL(num_eliminated, 1);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE_END()






