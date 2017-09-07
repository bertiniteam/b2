//This file is part of Bertini 2.
//
//pool_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//pool_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with pool_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire


//pool_test.cpp
//


 
#define BOOST_TEST_DYN_LINK 1

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 Object Pool Testing"
#include <boost/test/unit_test.hpp>


#include "logging.hpp"


using sec_level = boost::log::trivial::severity_level;

using LoggingInit = bertini::LoggingInit;


#include "bertini2/pool/system.hpp"


BOOST_GLOBAL_FIXTURE( LoggingInit );

BOOST_AUTO_TEST_SUITE(system_pool)

using namespace bertini;

BOOST_AUTO_TEST_CASE(make_system_pool)
{
	SystemPool sp;
}



BOOST_AUTO_TEST_CASE(make_nonpointer_system_and_add_to_pool)
{
	SystemPool sp;

	System sys;
	auto x = MakeVariable("x");
	auto y = MakeVariable("y");
	auto z = MakeVariable("z");
	
	sys.AddVariableGroup(VariableGroup({x,y,z}));  
	sys.AddFunction(x);
	sys.AddFunction(y);
	sys.AddFunction(z);

	auto result = sp.NonPtrAdd(sys);

	BOOST_CHECK(result.get() != &sys);
}


BOOST_AUTO_TEST_CASE(make_new_sys_from_pool)
{
	SystemPool sp;
	std::shared_ptr<System> sys = sp.NewObj();

	auto x = MakeVariable("x");
	auto y = MakeVariable("y");
	auto z = MakeVariable("z");
	
	sys->AddVariableGroup(VariableGroup({x,y,z}));  
	sys->AddFunction(x);
	sys->AddFunction(y);
	sys->AddFunction(z);
}


BOOST_AUTO_TEST_CASE(add_ptr_sys_to_pool)
{
	SystemPool sp;
	std::shared_ptr<System> sys = std::make_shared<System>();

	auto x = MakeVariable("x");
	auto y = MakeVariable("y");
	auto z = MakeVariable("z");
	
	sys->AddVariableGroup(VariableGroup({x,y,z}));  
	sys->AddFunction(x);
	sys->AddFunction(y);
	sys->AddFunction(z);

	auto result = sp.PtrAdd(sys);

	BOOST_CHECK(result.get() == sys.get());
}

BOOST_AUTO_TEST_SUITE_END()







BOOST_AUTO_TEST_SUITE(double_point_pool)

using namespace bertini;

BOOST_AUTO_TEST_CASE(make_double_point_point)
{
	PointPool<dbl> pool;
}
BOOST_AUTO_TEST_SUITE_END()




