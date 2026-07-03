// Parity + A/B benchmark for bertini::linalg::PartialPivLU (the stateful, allocation-lean
// multiprecision LU) against Eigen::PartialPivLU.  The parity suite runs in the normal test
// run; the timing suite is skipped unless BERTINI2_BENCH is set in the environment.

#include <boost/test/unit_test.hpp>

#include "bertini2/linalg/lu_solver.hpp"
#include "bertini2/mpfr_complex.hpp"

#include <Eigen/LU>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

using complex_dbl = bertini::complex_dbl;
using complex_mp  = bertini::complex_mp;
template<typename T> using Vec = bertini::Vec<T>;
template<typename T> using Mat = bertini::Mat<T>;

namespace {
	using Clock = std::chrono::steady_clock;

	// deterministic, well-conditioned (diagonally dominant) n x n matrix.
	template<typename NumT>
	Mat<NumT> MakeMatrix(unsigned n)
	{
		Mat<NumT> A(n, n);
		for (unsigned i = 0; i < n; ++i)
			for (unsigned j = 0; j < n; ++j)
			{
				double re = 0.3 + 0.017 * double((i * 7 + j * 3) % 11);
				double im = -0.2 + 0.011 * double((i * 5 + j * 13) % 9);
				A(i, j) = NumT(re, im);
			}
		for (unsigned i = 0; i < n; ++i)
			A(i, i) += NumT(double(n), 0.0);   // diagonal dominance -> well conditioned
		return A;
	}

	template<typename NumT>
	Vec<NumT> MakeRHS(unsigned n)
	{
		Vec<NumT> b(n);
		for (unsigned i = 0; i < n; ++i)
			b(i) = NumT(0.5 + 0.1 * double(i), 0.3 - 0.05 * double(i));
		return b;
	}
}

BOOST_AUTO_TEST_SUITE(LU_solver)

// Custom solver must agree with Eigen's PartialPivLU (small residual, matching solution).
BOOST_AUTO_TEST_CASE(parity_double)
{
	for (unsigned n : {3u, 5u, 20u, 50u})
	{
		auto A = MakeMatrix<complex_dbl>(n);
		auto b = MakeRHS<complex_dbl>(n);

		Eigen::PartialPivLU<Mat<complex_dbl>> elu(A);
		Vec<complex_dbl> xe = elu.solve(b);

		bertini::linalg::PartialPivLU<complex_dbl> blu;
		blu.ChangeSize(n);
		BOOST_CHECK(blu.Factor(A) == bertini::MatrixSuccessCode::Success);
		Vec<complex_dbl> xc(n);
		blu.Solve(b, xc);

		double resid = (A * xc - b).norm();
		double diff  = (xc - xe).norm();
		BOOST_CHECK_SMALL(resid, 1e-12);
		BOOST_CHECK_SMALL(diff,  1e-12);
	}
}

BOOST_AUTO_TEST_CASE(parity_mp)
{
	auto saved = bertini::DefaultPrecision();
	bertini::DefaultPrecision(40);

	for (unsigned n : {3u, 5u, 20u, 50u})
	{
		auto A = MakeMatrix<complex_mp>(n);
		auto b = MakeRHS<complex_mp>(n);

		Eigen::PartialPivLU<Mat<complex_mp>> elu(A);
		Vec<complex_mp> xe = elu.solve(b);

		bertini::linalg::PartialPivLU<complex_mp> blu;
		blu.ChangeSize(n);
		blu.ChangePrecision(40);
		BOOST_CHECK(blu.Factor(A) == bertini::MatrixSuccessCode::Success);
		Vec<complex_mp> xc(n);
		blu.Solve(b, xc);

		double resid = static_cast<double>((A * xc - b).norm());
		double diff  = static_cast<double>((xc - xe).norm());
		BOOST_CHECK_SMALL(resid, 1e-30);
		BOOST_CHECK_SMALL(diff,  1e-30);
	}

	bertini::DefaultPrecision(saved);
}

// A/B timing: Eigen vs custom, factor+solve, ns/op.  Skipped unless BERTINI2_BENCH is set.
BOOST_AUTO_TEST_CASE(ab_timing)
{
	if (!std::getenv("BERTINI2_BENCH")) return;

	bertini::DefaultPrecision(40);
	std::cout << "\nLU_BENCH_BEGIN  (factor+solve, complex_mp @ 40 digits; ns/op, lower=better)\n";
	std::cout << "| N | eigen ns | custom ns | speedup |\n|---|---:|---:|---:|\n";

	struct Case { unsigned n; std::size_t iters; };
	for (Case c : std::vector<Case>{ {5, 20000}, {20, 2000}, {50, 300}, {100, 60} })
	{
		unsigned n = c.n;
		auto A = MakeMatrix<complex_mp>(n);
		auto b = MakeRHS<complex_mp>(n);
		Vec<complex_mp> x(n);

		// Eigen
		Eigen::PartialPivLU<Mat<complex_mp>> elu(n);
		elu.compute(A); x = elu.solve(b);                       // warmup
		auto t0 = Clock::now();
		for (std::size_t m = 0; m < c.iters; ++m) { elu.compute(A); x = elu.solve(b); }
		auto t1 = Clock::now();
		double eigen_ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / double(c.iters);

		// custom
		bertini::linalg::PartialPivLU<complex_mp> blu;
		blu.ChangeSize(n); blu.ChangePrecision(40);
		blu.Factor(A); blu.Solve(b, x);                          // warmup
		t0 = Clock::now();
		for (std::size_t m = 0; m < c.iters; ++m) { blu.Factor(A); blu.Solve(b, x); }
		t1 = Clock::now();
		double custom_ns = std::chrono::duration<double, std::nano>(t1 - t0).count() / double(c.iters);

		std::cout << "| " << n
		          << " | " << std::fixed << std::setprecision(0) << eigen_ns
		          << " | " << custom_ns
		          << " | " << std::setprecision(2) << (eigen_ns / custom_ns) << "x |\n";
	}
	std::cout << "LU_BENCH_END\n";
}

BOOST_AUTO_TEST_SUITE_END()
