//This file is part of Bertini 2.
//
//degrees_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//degrees_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this file.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file degrees_test.cpp

\brief Degrees of block-composed Systems.

Now that the polynomial path is folded into a PolynomialBlock and every structured block
(products-of-linears, linear-forms, blend) reports its own degrees, System::Degrees() /
Degrees(vars) / DegreeBound() concatenate degrees across all blocks.  These tests pin the
per-block degree semantics -- with particular attention to the linear-algebra path
(LinearFormsBlock = degree 1), where being seen as degree-1 is exactly what makes A@x+b a
linear form to the start-system / degree machinery.

Block degree contract under test:
  PolynomialBlock      -- the function tree's Degree() / Degree(vars)
  LinearFormsBlock     -- always 1
  ProductsOfLinearsBlock -- the number of linear factors (== rows of the factor matrix)
  BlendBlock           -- the elementwise max of its operands' degrees
*/

#include <boost/test/unit_test.hpp>

#include "bertini2/system/system.hpp"
#include "bertini2/system/start_systems.hpp"
#include "bertini2/system/blocks/block.hpp"
#include "bertini2/system/blocks/polynomial_block.hpp"
#include "bertini2/system/blocks/products_of_linears_block.hpp"
#include "bertini2/system/blocks/linear_forms_block.hpp"
#include "bertini2/system/blocks/blend_block.hpp"

BOOST_AUTO_TEST_SUITE(degrees_suite)

using namespace bertini;
using Var = std::shared_ptr<node::Variable>;
using bertini::node::Variable;
using bertini::node::Integer;
using bertini::blocks::PolynomialBlock;
using bertini::blocks::ProductsOfLinearsBlock;
using bertini::blocks::LinearFormsBlock;
using bertini::blocks::BlendBlock;
using bertini::DefaultPrecision;

namespace {

// augmented coefficient matrix (rows x (num_vars+1)), last column the constant term.
Mat<complex_mp> AugMat(std::vector<std::vector<int>> const& rows)
{
	Mat<complex_mp> M(static_cast<Eigen::Index>(rows.size()),
	                    static_cast<Eigen::Index>(rows.front().size()));
	for (Eigen::Index i = 0; i < M.rows(); ++i)
		for (Eigen::Index j = 0; j < M.cols(); ++j)
			M(i, j) = complex_mp(rows[static_cast<size_t>(i)][static_cast<size_t>(j)]);
	return M;
}

} // namespace


////////////////////////////////////////////////////////////////////////////////
// LinearFormsBlock -- the linear-algebra path.  Every function is degree 1.
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(linear_forms_block_is_degree_one)
{
	DefaultPrecision(30);
	// three linear forms in two variables (A @ x + b shape): 3 x (2+1) augmented matrix.
	LinearFormsBlock block(2, AugMat({{2, 3, 1}, {1, -1, 4}, {0, 5, -2}}));

	std::vector<int> expected{1, 1, 1};
	BOOST_CHECK(block.Degrees() == expected);

	// the structured-block contract: Degrees(group) forwards to the total Degrees().
	VariableGroup grp{Variable::Make("x"), Variable::Make("y")};
	BOOST_CHECK(block.Degrees(grp) == block.Degrees());
}

BOOST_AUTO_TEST_CASE(system_of_only_linear_forms_has_degree_bound_one)
{
	DefaultPrecision(30);
	System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddBlock(LinearFormsBlock(2, AugMat({{2, 3, 1}, {1, -1, 4}})));

	std::vector<int> expected{1, 1};
	BOOST_CHECK(sys.Degrees() == expected);
	BOOST_CHECK(sys.Degrees(sys.Variables()) == expected);
	BOOST_CHECK_EQUAL(sys.DegreeBound(), 1);
}


////////////////////////////////////////////////////////////////////////////////
// PolynomialBlock -- the folded classic path.  Degrees come from the trees.
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(polynomial_block_degrees_from_trees)
{
	DefaultPrecision(30);
	System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x * y - Integer::Make(1));     // degree 2
	sys.AddFunction(x + y);                         // degree 1
	sys.AddFunction(Integer::Make(7));              // degree 0 (constant)

	std::vector<int> expected{2, 1, 0};
	BOOST_CHECK(sys.Degrees() == expected);
	BOOST_CHECK_EQUAL(sys.DegreeBound(), 2);
}

BOOST_AUTO_TEST_CASE(polynomial_degrees_match_hand_derived)
{
	// behavior-preservation guard: folding functions into the PolynomialBlock did not change
	// degree reporting -- it is still the function tree's Degree().
	DefaultPrecision(30);
	System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x * x * y * y * y);             // x^2 y^3, total degree 5
	sys.AddFunction(x + y + Integer::Make(1));      // degree 1

	std::vector<int> expected{5, 1};
	BOOST_CHECK(sys.Degrees(sys.Variables()) == expected);
	BOOST_CHECK_EQUAL(sys.DegreeBound(), 5);
}


////////////////////////////////////////////////////////////////////////////////
// The eigenvalue / linear-algebra distinction: (A - lam I) x is *bilinear* (degree 2),
// while a constant-coefficient normalization c.x - 1 is *linear* (degree 1).  The degree
// machinery must see exactly this difference.
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(bilinear_eigenvalue_row_is_degree_two)
{
	DefaultPrecision(30);
	Var x0 = Variable::Make("x0"), x1 = Variable::Make("x1"), lam = Variable::Make("lam");
	// row of (A - lam I) x : a*x0 + b*x1 - lam*x0   (the lam*x0 term is degree 2)
	auto row = Integer::Make(2) * x0 + Integer::Make(3) * x1 - lam * x0;
	BOOST_CHECK_EQUAL(row->Degree(), 2);
}

BOOST_AUTO_TEST_CASE(eigenproblem_degrees_bilinear_rows_plus_linear_normalization)
{
	DefaultPrecision(30);
	System sys;
	Var x0 = Variable::Make("x0"), x1 = Variable::Make("x1"), lam = Variable::Make("lam");
	sys.AddVariableGroup(VariableGroup{x0, x1, lam});

	// the two (A - lam I) x rows -- degree 2 each (eigenvalue * eigenvector entry).
	sys.AddFunction(Integer::Make(4) * x0 + Integer::Make(1) * x1 - lam * x0);
	sys.AddFunction(Integer::Make(1) * x0 + Integer::Make(3) * x1 - lam * x1);
	// generic normalization c.x - 1 as a linear-forms block -- degree 1.
	// (LinearFormsBlock spans all 3 columns; the lam coefficient is 0.)
	sys.AddBlock(LinearFormsBlock(3, AugMat({{2, 5, 0, -1}})));

	// poly block (added first by AddFunction) then the linear-forms block.
	std::vector<int> expected{2, 2, 1};
	BOOST_CHECK(sys.Degrees() == expected);
	BOOST_CHECK(sys.Degrees(sys.Variables()) == expected);
	BOOST_CHECK_EQUAL(sys.DegreeBound(), 2);
}


////////////////////////////////////////////////////////////////////////////////
// ProductsOfLinearsBlock -- degree == number of linear factors.
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(products_of_linears_degree_is_factor_count)
{
	DefaultPrecision(30);
	// f0 = product of 2 linear factors (degree 2); f1 = product of 3 (degree 3).
	Mat<complex_mp> f0 = AugMat({{2, 3, 1}, {1, -1, 4}});
	Mat<complex_mp> f1 = AugMat({{1, 0, 1}, {0, 1, -2}, {1, 1, 0}});
	ProductsOfLinearsBlock block(2, std::vector<Mat<complex_mp>>{f0, f1});

	std::vector<int> expected{2, 3};
	BOOST_CHECK(block.Degrees() == expected);
	// structured-block contract: Degrees(group) forwards to total Degrees().
	VariableGroup grp{Variable::Make("x"), Variable::Make("y")};
	BOOST_CHECK(block.Degrees(grp) == block.Degrees());

	System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddBlock(ProductsOfLinearsBlock(2, std::vector<Mat<complex_mp>>{f0, f1}));
	BOOST_CHECK(sys.Degrees() == expected);
	BOOST_CHECK_EQUAL(sys.DegreeBound(), 3);
}


////////////////////////////////////////////////////////////////////////////////
// BlendBlock -- elementwise max of its operands' degrees (the homotopy case).
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(blend_degree_is_elementwise_max_of_operands)
{
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");

	auto target = std::make_shared<System>();
	target->AddVariableGroup(VariableGroup{x, y});
	target->AddFunction(x * y);                    // degree 2
	target->AddFunction(x * x * x);                // degree 3

	auto start = std::make_shared<System>();
	start->AddVariableGroup(VariableGroup{x, y});
	start->AddFunction(x);                          // degree 1
	start->AddFunction(y);                          // degree 1

	std::vector<std::shared_ptr<node::Node>> coeffs{ Integer::Make(1) - t, t };
	std::vector<std::shared_ptr<const System>> operands{ target, start };
	BlendBlock<System> blend(t, std::move(coeffs), operands);

	std::vector<int> expected{2, 3};               // max({2,3},{1,1})
	BOOST_CHECK(blend.Degrees() == expected);
}

BOOST_AUTO_TEST_CASE(system_with_blend_block_has_nonempty_degree_bound)
{
	// regression: the MHom-style blend homotopy carries no polynomial block, yet its degrees
	// (and DegreeBound, needed by the AMP config) must come from the blend's operands rather
	// than being empty (which previously dereferenced an empty range in DegreeBound()).
	DefaultPrecision(30);
	Var x = Variable::Make("x"), y = Variable::Make("y"), t = Variable::Make("t");

	auto target = std::make_shared<System>();
	target->AddVariableGroup(VariableGroup{x, y});
	target->AddFunction(x * y);                    // degree 2
	target->AddFunction(x * x * x);                // degree 3

	auto start = std::make_shared<System>();
	start->AddVariableGroup(VariableGroup{x, y});
	start->AddFunction(x);
	start->AddFunction(y);

	System H;
	H.AddVariableGroup(VariableGroup{x, y});
	H.AddPathVariable(t);
	std::vector<std::shared_ptr<node::Node>> coeffs{ Integer::Make(1) - t, t };
	std::vector<std::shared_ptr<const System>> operands{ target, start };
	H.AddBlock(BlendBlock<System>(t, std::move(coeffs), operands));

	BOOST_REQUIRE(H.HasStructuredBlocks());             // it does have a structured block
	std::vector<int> expected{2, 3};
	BOOST_CHECK(H.Degrees() == expected);
	BOOST_CHECK_EQUAL(H.DegreeBound(), 3);
}


////////////////////////////////////////////////////////////////////////////////
// Mixed and edge cases.
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(mixed_poly_and_linear_forms_concatenate_in_block_order)
{
	DefaultPrecision(30);
	System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x, y});
	sys.AddFunction(x * x);                         // poly block, degree 2
	sys.AddFunction(x * x * x);                     // poly block, degree 3
	sys.AddBlock(LinearFormsBlock(2, AugMat({{2, 3, 1}, {1, -1, 4}})));  // degree 1, 1

	// AddFunction creates the polynomial block first, then the linear-forms block is appended.
	std::vector<int> expected{2, 3, 1, 1};
	BOOST_CHECK(sys.Degrees() == expected);
	BOOST_CHECK_EQUAL(sys.DegreeBound(), 3);
}

BOOST_AUTO_TEST_CASE(system_with_no_functions_has_no_degrees)
{
	// A System with variables but no functions/blocks reports no degrees.  (DegreeBound() has
	// a non-empty precondition -- it max-reduces the degree list -- so it is not called here.)
	DefaultPrecision(30);
	System sys;
	sys.AddVariableGroup(VariableGroup{Variable::Make("x"), Variable::Make("y")});
	BOOST_CHECK(sys.Degrees().empty());
	BOOST_CHECK(sys.Degrees(sys.Variables()).empty());
}


////////////////////////////////////////////////////////////////////////////////
// The real MHom start system: products-of-linears block, degrees == factor counts,
// and (the regression) a finite DegreeBound.
////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(mhom_start_system_has_factor_count_degrees)
{
	DefaultPrecision(30);
	using namespace bertini::start_system;

	System sys;
	Var x = Variable::Make("x"), y = Variable::Make("y");
	sys.AddVariableGroup(VariableGroup{x});
	sys.AddVariableGroup(VariableGroup{y});
	sys.AddFunction(x * y - Integer::Make(1));
	sys.AddFunction(x + y);
	sys.Homogenize();
	sys.AutoPatch();

	auto mhom = MHomogeneous(sys);

	auto degs = mhom.Degrees();
	BOOST_REQUIRE(!degs.empty());                                  // was empty -> crashed DegreeBound
	BOOST_CHECK_EQUAL(degs.size(), mhom.NumNaturalFunctions());
	for (int d : degs)
		BOOST_CHECK_GE(d, 1);                                      // each start function is a product of >=1 linear
	BOOST_CHECK_GE(mhom.DegreeBound(), 1);                         // finite, no empty-range deref
}

BOOST_AUTO_TEST_SUITE_END()
