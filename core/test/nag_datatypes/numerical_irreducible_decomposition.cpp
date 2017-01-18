//This file is part of Bertini 2.
//
//test/nag_algorithms/numerical_irreducible_decomposition.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//test/nag_algorithms/numerical_irreducible_decomposition.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with test/nag_algorithms/numerical_irreducible_decomposition.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

/**
\file test/nag_datatypes/numerical_irreducible_decomposition.cpp  Tests the NID algorithm.
*/

// individual authors of this file include:
// dani brake, university of notre dame

#include <boost/test/unit_test.hpp>
#include "bertini2/nag_datatypes/numerical_irreducible_decomposition.hpp"
#include "bertini2/system/precon.hpp"



BOOST_AUTO_TEST_SUITE(nid)



	BOOST_AUTO_TEST_SUITE(policy_by_shared_pointer)

		template <typename T>
		using sp = std::shared_ptr<T>;

		using NID = bertini::nag_datatype::NumericalIrreducibleDecomposition<
			bertini::complex>;

		BOOST_AUTO_TEST_CASE(something)
		{	

		}

	BOOST_AUTO_TEST_SUITE_END() // by_shared_pointer


BOOST_AUTO_TEST_SUITE_END() // witness_set

