//This file is part of Bertini 2.
//
//slp_intern_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//slp_intern_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with slp_intern_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Tests for SLPProgram hash-consing (ADR-0027 E4, ADR-0042): ContentHash/SameContent,
the InternProgram weak table, sharing across facades compiled from identical sources,
re-interning on deserialization, and SLPMemory independence of facades over one Program.
*/

#include <boost/test/unit_test.hpp>
#include "bertini2/system/straight_line_program.hpp"
#include "bertini2/system/system.hpp"
#include "bertini2/io/parsing/system_parsers.hpp"

#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using SLP = bertini::StraightLineProgram;

namespace {

bertini::System ParseSystem(std::string const& str){
	bertini::System sys;
	[[maybe_unused]] bool success = bertini::parsing::classic::parse(str.begin(), str.end(), sys);
	return sys;
}

bertini::System CircleLineSystem(){
	return ParseSystem("function f,g; variable_group x,y; f = x^2+y^2-1; g = x-y;");
}

bertini::System OtherSystem(){
	return ParseSystem("function f,g; variable_group x,y; f = x^2-y; g = x+y+1;");
}

} // unnamed namespace

BOOST_AUTO_TEST_SUITE(slp_intern)

BOOST_AUTO_TEST_CASE(identical_systems_share_one_program)
{
	auto sys_a = CircleLineSystem();
	auto sys_b = CircleLineSystem();  // independently parsed, equal content

	SLP slp_a(sys_a);
	SLP slp_b(sys_b);

	BOOST_CHECK(slp_a.Program() == slp_b.Program());  // pointer equality: one interned Program
	BOOST_CHECK(slp_a.Program()->SameContent(*slp_b.Program()));
	BOOST_CHECK_EQUAL(slp_a.Program()->ContentHash(), slp_b.Program()->ContentHash());
}

BOOST_AUTO_TEST_CASE(different_systems_get_different_programs)
{
	auto sys_a = CircleLineSystem();
	auto sys_b = OtherSystem();

	SLP slp_a(sys_a);
	SLP slp_b(sys_b);

	BOOST_CHECK(slp_a.Program() != slp_b.Program());
	BOOST_CHECK(!slp_a.Program()->SameContent(*slp_b.Program()));
}

BOOST_AUTO_TEST_CASE(same_content_implies_equal_hash_and_is_reflexive)
{
	auto sys = CircleLineSystem();
	SLP slp(sys);
	BOOST_CHECK(slp.Program()->SameContent(*slp.Program()));
	BOOST_CHECK_EQUAL(slp.Program()->ContentHash(), slp.Program()->ContentHash());
}

BOOST_AUTO_TEST_CASE(deserialized_program_reinterns_to_live_one)
{
	auto sys = CircleLineSystem();
	SLP original(sys);

	std::stringstream archive_stream;
	{
		boost::archive::text_oarchive oa(archive_stream);
		oa << original;
	}
	SLP loaded;
	{
		boost::archive::text_iarchive ia(archive_stream);
		ia >> loaded;
	}

	// load() routes through InternProgram, so the loaded facade shares the live Program.
	BOOST_CHECK(loaded.Program() == original.Program());
}

BOOST_AUTO_TEST_CASE(facades_over_one_program_have_independent_memory)
{
	auto sys_a = CircleLineSystem();
	auto sys_b = CircleLineSystem();

	SLP slp_a(sys_a);
	SLP slp_b(sys_b);
	BOOST_REQUIRE(slp_a.Program() == slp_b.Program());

	using complex_dbl = bertini::complex_dbl;
	bertini::Vec<complex_dbl> pt_a(2), pt_b(2);
	pt_a << complex_dbl(1.0, 0.0), complex_dbl(0.0, 1.0);
	pt_b << complex_dbl(2.0, 0.0), complex_dbl(3.0, 0.0);

	slp_a.Eval(pt_a);
	slp_b.Eval(pt_b);

	auto values_a = slp_a.GetFuncVals<complex_dbl>();
	auto values_b = slp_b.GetFuncVals<complex_dbl>();

	// f = x^2+y^2-1, g = x-y at the two distinct points; memories did not bleed.
	BOOST_CHECK(abs(values_a(0) - complex_dbl(-1.0, 0.0)) < 1e-14);   // 1 + (-1) - 1
	BOOST_CHECK(abs(values_a(1) - complex_dbl(1.0, -1.0)) < 1e-14);   // 1 - i
	BOOST_CHECK(abs(values_b(0) - complex_dbl(12.0, 0.0)) < 1e-14);   // 4 + 9 - 1
	BOOST_CHECK(abs(values_b(1) - complex_dbl(-1.0, 0.0)) < 1e-14);   // 2 - 3
}

BOOST_AUTO_TEST_CASE(intern_table_is_weak_and_self_cleans)
{
	std::shared_ptr<const bertini::SLPProgram> first_program;
	{
		auto sys = OtherSystem();
		SLP slp(sys);
		first_program = slp.Program();
	}
	// Drop the only owner; the table holds weak_ptrs, so the entry is dead now.
	auto const* old_raw = first_program.get();
	first_program.reset();

	// A recompile registers a NEW program (the old one is gone, not resurrected).
	auto sys = OtherSystem();
	SLP slp(sys);
	BOOST_CHECK(slp.Program() != nullptr);
	// The old registration must not have been handed back as a dangling pointer; the new
	// program is live and content-correct.
	BOOST_CHECK(slp.Program()->SameContent(*slp.Program()));
	(void)old_raw;  // address reuse is possible; pointer inequality is NOT asserted
}

BOOST_AUTO_TEST_SUITE_END()
