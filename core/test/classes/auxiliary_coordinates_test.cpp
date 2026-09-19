//This file is part of Bertini 2.
//
//auxiliary_coordinates_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//auxiliary_coordinates_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with auxiliary_coordinates_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file Auxiliary coordinates: the ones a system carries for the construction's sake, which every
judgement about a point leaves out (b2#403).

A system in null-vector form -- `f; v^T M; patch` -- carries `v` with its own patch, necessarily,
or the system would admit the zero null vector.  That patch fixes a normalization nobody chose on
geometric grounds, so `v`'s DIRECTION is the meaningful object and its MAGNITUDE is an artifact:
not intentionally divergent, simply not governed.  Judging "is this point finite" on the larger of
two independently uncontrolled scales is a kind error rather than a tuning problem -- there is no
correct threshold for an ungoverned quantity.

So the system is told which of its coordinates are auxiliary, and it answers the questions itself.
Bertini attaches no further meaning to the word: an auxiliary coordinate is still tracked, recorded
and returned, it is simply not evidence about whether the point ran off to infinity.
*/

#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "bertini2/system/system.hpp"
#include "bertini2/system/start_systems.hpp"

using namespace bertini;
using Variable = node::Variable;

namespace {

// (x) and (v): two affine groups, the second standing in for a block whose scale nobody chose.
System TwoGroups()
{
    auto x = Variable::Make("x");
    auto v = Variable::Make("v");
    System sys;
    sys.AddVariableGroup(VariableGroup{x});
    sys.AddVariableGroup(VariableGroup{v});
    sys.AddFunction(pow(x, 2) - 1);
    sys.AddFunction(v - 1);
    return sys;
}

Vec<complex_dbl> Point(complex_dbl a, complex_dbl b)
{
    Vec<complex_dbl> p(2);
    p(0) = a;
    p(1) = b;
    return p;
}

} // unnamed namespace


BOOST_AUTO_TEST_SUITE(auxiliary_coordinates)

BOOST_AUTO_TEST_CASE(by_default_every_coordinate_is_judged)
{
    auto sys = TwoGroups();
    BOOST_CHECK(!sys.HaveAuxiliaryCoordinates());
    BOOST_CHECK_EQUAL(sys.NumAuxiliaryCoordinates(), 0u);

    BOOST_CHECK(sys.IsFinite(Point(1.0, 1.0), 1e5));
    BOOST_CHECK(!sys.IsFinite(Point(1.0, 1e9), 1e5));   // the second coordinate carries the verdict
    BOOST_CHECK(!sys.IsFinite(Point(1e9, 1.0), 1e5));
}


BOOST_AUTO_TEST_CASE(an_auxiliary_group_is_not_evidence_about_infinity)
{
    auto sys = TwoGroups();
    sys.SetAuxiliaryVariableGroups({1});

    BOOST_CHECK(sys.HaveAuxiliaryCoordinates());
    BOOST_CHECK_EQUAL(sys.NumAuxiliaryCoordinates(), 1u);

    // the case from the issue: natural coordinates well inside the working region, the ungoverned
    // block far outside it, and the point is finite because the block was never the question
    BOOST_CHECK(sys.IsFinite(Point(1e3, 1e9), 1e5));
    // and it excludes rather than disables: the coordinate still being judged decides
    BOOST_CHECK(!sys.IsFinite(Point(1e9, 1.0), 1e5));
}


BOOST_AUTO_TEST_CASE(individual_coordinates_can_be_auxiliary_when_the_grouping_does_not_separate_them)
{
    auto x = Variable::Make("x");
    auto v = Variable::Make("v");
    System sys;
    sys.AddVariableGroup(VariableGroup{x, v});    // deliberately one group
    sys.AddFunction(pow(x, 2) - 1);
    sys.AddFunction(v - 1);

    sys.SetAuxiliaryCoordinates({1});
    BOOST_CHECK_EQUAL(sys.NumAuxiliaryCoordinates(), 1u);
    BOOST_CHECK(sys.IsFinite(Point(1e3, 1e9), 1e5));
    BOOST_CHECK(!sys.IsFinite(Point(1e9, 1.0), 1e5));
}


BOOST_AUTO_TEST_CASE(auxiliary_groups_and_auxiliary_coordinates_are_a_union)
{
    auto x = Variable::Make("x");
    auto y = Variable::Make("y");
    auto v = Variable::Make("v");
    System sys;
    sys.AddVariableGroup(VariableGroup{x, y});
    sys.AddVariableGroup(VariableGroup{v});
    sys.AddFunction(x + y + v);

    sys.SetAuxiliaryVariableGroups({1});      // v
    sys.SetAuxiliaryCoordinates({1});         // y
    BOOST_CHECK_EQUAL(sys.NumAuxiliaryCoordinates(), 2u);
    BOOST_CHECK(!sys.CoordinateIsAuxiliary(0));
    BOOST_CHECK(sys.CoordinateIsAuxiliary(1));
    BOOST_CHECK(sys.CoordinateIsAuxiliary(2));
}


BOOST_AUTO_TEST_CASE(realness_is_judged_over_the_same_coordinates_as_finiteness)
{
    // the identical defect, one line away in the classifier: a point whose only complex
    // coordinate is one the system has been told to disregard is real
    auto sys = TwoGroups();
    BOOST_CHECK(!sys.IsReal(Point(1.0, complex_dbl(1.0, 0.5)), 1e-8));

    sys.SetAuxiliaryVariableGroups({1});
    BOOST_CHECK(sys.IsReal(Point(1.0, complex_dbl(1.0, 0.5)), 1e-8));
    BOOST_CHECK(!sys.IsReal(Point(complex_dbl(1.0, 0.5), 1.0), 1e-8));
}


BOOST_AUTO_TEST_CASE(auxiliary_coordinates_survive_homogenization_and_patching)
{
    // the indices are into a point in USER coordinates, so adding a homogenizing variable per
    // group and rescaling onto a patch does not move them -- which is what lets a declaration
    // made on the system the author wrote still mean the same thing inside the solver
    auto sys = TwoGroups();
    sys.SetAuxiliaryVariableGroups({1});
    sys.Homogenize();
    sys.AutoPatch();

    BOOST_REQUIRE(sys.IsPatched());
    auto const internal = sys.HomogenizePoint(Point(1e3, 1e9));
    BOOST_CHECK(sys.IsFinite(internal, 1e5));
    BOOST_CHECK(!sys.IsFinite(sys.HomogenizePoint(Point(1e9, 1.0)), 1e5));
}


BOOST_AUTO_TEST_CASE(a_copy_is_judged_the_way_its_original_is)
{
    // System's copy constructor is hand-written, member by member -- a mirror like serialize and
    // the canonical encoder -- so a member it forgets is silently dropped by every copy.  It
    // forgot these, and since Clone IS a copy, the solver's target arrived with an empty set and
    // the whole feature did nothing past the system the caller was holding.
    auto sys = TwoGroups();
    sys.SetAuxiliaryVariableGroups({1});

    auto const copied = sys;
    BOOST_CHECK(copied.AuxiliaryVariableGroups() == sys.AuxiliaryVariableGroups());
    BOOST_CHECK(Clone(sys).AuxiliaryVariableGroups() == sys.AuxiliaryVariableGroups());
    BOOST_CHECK(Clone(sys).IsFinite(Point(1e3, 1e9), 1e5));
}


BOOST_AUTO_TEST_CASE(a_homotopy_is_judged_the_way_its_target_is)
{
    // the tracker measures the HOMOTOPY, not the system the author wrote, so auxiliary
    // coordinates that did not travel would be silently ignored for the whole of tracking
    auto target = TwoGroups();
    target.SetAuxiliaryVariableGroups({1});

    auto start = TwoGroups();              // same variable structure, different equations
    auto homotopy = MakeHomotopy(target, start, "t");

    BOOST_CHECK(homotopy.AuxiliaryVariableGroups() == target.AuxiliaryVariableGroups());
    BOOST_CHECK(homotopy.AuxiliaryCoordinates() == target.AuxiliaryCoordinates());
}


BOOST_AUTO_TEST_CASE(which_coordinates_are_auxiliary_is_part_of_what_the_system_is)
{
    // it changes what gets tracked -- the tracker truncates on what is NOT auxiliary -- so two
    // systems differing only here are asking different questions and must not recall each other
    auto plain = TwoGroups();
    auto with_aux = TwoGroups();
    with_aux.SetAuxiliaryVariableGroups({1});

    BOOST_CHECK(plain.ContentDigest() != with_aux.ContentDigest());

    auto same = TwoGroups();
    same.SetAuxiliaryVariableGroups({1});
    BOOST_CHECK(same.ContentDigest() == with_aux.ContentDigest());
}


BOOST_AUTO_TEST_CASE(the_order_they_were_given_in_is_not_part_of_the_identity)
{
    auto x = Variable::Make("x");
    auto y = Variable::Make("y");
    auto z = Variable::Make("z");
    auto build = [&](std::vector<unsigned> auxiliary) {
        System sys;
        sys.AddVariableGroup(VariableGroup{x, y, z});
        sys.AddFunction(x + y + z);
        sys.SetAuxiliaryCoordinates(std::move(auxiliary));
        return sys;
    };

    BOOST_CHECK(build({1, 2}).ContentDigest() == build({2, 1}).ContentDigest());
    BOOST_CHECK(build({1, 2}).ContentDigest() == build({2, 1, 2}).ContentDigest());
}


BOOST_AUTO_TEST_CASE(auxiliary_coordinates_round_trip_through_the_canonical_encoding)
{
    auto x = Variable::Make("x");
    auto y = Variable::Make("y");
    auto v = Variable::Make("v");
    System sys;
    sys.AddVariableGroup(VariableGroup{x, y});
    sys.AddVariableGroup(VariableGroup{v});
    sys.AddFunction(x + y + v);
    sys.SetAuxiliaryVariableGroups({1});      // v
    sys.SetAuxiliaryCoordinates({1});         // y

    auto const rebuilt = System::FromCanonicalEncoding(sys.CanonicalEncodingText());
    BOOST_CHECK(rebuilt.AuxiliaryVariableGroups() == sys.AuxiliaryVariableGroups());
    BOOST_CHECK(rebuilt.AuxiliaryCoordinates() == sys.AuxiliaryCoordinates());
    BOOST_CHECK(rebuilt.ContentDigest() == sys.ContentDigest());
}


BOOST_AUTO_TEST_CASE(an_index_that_names_nothing_is_refused)
{
    auto sys = TwoGroups();
    BOOST_CHECK_THROW(sys.SetAuxiliaryVariableGroups({2}), std::out_of_range);
    BOOST_CHECK_THROW(sys.SetAuxiliaryCoordinates({2}), std::out_of_range);
    BOOST_CHECK(!sys.HaveAuxiliaryCoordinates());   // a refused call changes nothing
}


BOOST_AUTO_TEST_CASE(making_every_coordinate_auxiliary_is_refused)
{
    // a point with nothing to judge would be unconditionally finite and real, which is not a
    // verdict this library should render
    auto sys = TwoGroups();
    BOOST_CHECK_THROW(sys.SetAuxiliaryVariableGroups({0, 1}), std::runtime_error);
    BOOST_CHECK_THROW(sys.SetAuxiliaryCoordinates({0, 1}), std::runtime_error);

    // and the refusal sees the UNION, so it cannot be walked around one setter at a time
    sys.SetAuxiliaryVariableGroups({1});
    BOOST_CHECK_THROW(sys.SetAuxiliaryCoordinates({0}), std::runtime_error);
    BOOST_CHECK_EQUAL(sys.NumAuxiliaryCoordinates(), 1u);   // the earlier declaration survives
}

BOOST_AUTO_TEST_SUITE_END()
