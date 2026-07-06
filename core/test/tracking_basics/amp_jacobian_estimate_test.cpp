//This file is part of Bertini 2.
//
//amp_jacobian_estimate_test.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//amp_jacobian_estimate_test.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with amp_jacobian_estimate_test.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file amp_jacobian_estimate_test.cpp

\brief Verifies the ||J^{-1}|| estimate the AMP corrector uses to drive precision.

The corrector estimates ||J^{-1}|| as ||J^{-1} r|| for a random unit-modulus vector r (one LU
solve), and feeds it into CriterionB -> DigitsB.  The AMP-escalation investigation found DigitsB
spiking to ~147 (forcing precision to ~140 digits).  Since DigitsB ~ log10(||J^{-1}|| * (...)),
that requires ||J^{-1}|| ~ 1e145 -- IF the estimate is faithful.  These tests check that the
estimate faithfully tracks the true spectral ||J^{-1}|| = 1/sigma_min(J): it is a lower bound up
to the ||r||=sqrt(n) factor, and it does NOT spuriously inflate by many orders of magnitude.  So
a reported huge ||J^{-1}|| reflects a genuine near-singularity, not an estimation artifact.
*/

#include <boost/test/unit_test.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include "bertini2/mpfr_complex.hpp"
#include "bertini2/eigen_extensions.hpp"

BOOST_AUTO_TEST_SUITE(amp_jacobian_estimate)

using complex_dbl = bertini::complex_dbl;
template <typename T> using Vec = bertini::Vec<T>;
template <typename T> using Mat = bertini::Mat<T>;

namespace {

	// The corrector's estimate: ||J^{-1} r|| for a random unit-modulus r (max over a few draws,
	// as a path accumulates many draws).  This mirrors newton_corrector.hpp.
	double NormJInverseEstimate(Mat<complex_dbl> const& J, int draws = 12)
	{
		double best = 0.0;
		auto lu = J.partialPivLu();
		for (int i = 0; i < draws; ++i)
		{
			Vec<complex_dbl> r = bertini::RandomOfUnits<complex_dbl>(static_cast<unsigned>(J.cols()));
			best = std::max(best, lu.solve(r).norm());
		}
		return best;
	}

	// True spectral ||J^{-1}|| = 1 / smallest singular value of J.
	double TrueNormJInverse(Mat<complex_dbl> const& J)
	{
		Eigen::JacobiSVD<Mat<complex_dbl>> svd(J);
		double sigma_min = svd.singularValues()(svd.singularValues().size() - 1);
		return 1.0 / sigma_min;
	}

} // anonymous namespace


// For well-conditioned matrices the estimate is the right order of magnitude (it never reports a
// huge ||J^{-1}|| for a benign matrix -- the spike is not manufactured here).
BOOST_AUTO_TEST_CASE(estimate_tracks_truth_well_conditioned)
{
	for (int n = 2; n <= 6; ++n)
		for (int trial = 0; trial < 20; ++trial)
		{
			Mat<complex_dbl> J = Mat<complex_dbl>::Random(n, n);
			double est  = NormJInverseEstimate(J);
			double tru  = TrueNormJInverse(J);

			// estimate is a lower bound up to the ||r|| = sqrt(n) factor: est <= sqrt(n)*true.
			BOOST_CHECK_LE(est, std::sqrt(double(n)) * tru * (1 + 1e-9));
			// and it captures the magnitude: not absurdly small either.
			BOOST_CHECK_GE(est, tru / (100.0 * n));
		}
}


// A deliberately near-singular matrix with sigma_min ~ 1e-k.  Representing this conditioning at
// all requires precision > k digits -- exactly why AMP escalates -- so this is done in mpfr.  The
// estimate must report ~1e+k: it faithfully reflects the genuine ill-conditioning, neither missing
// it nor inflating it by orders.  This is the case the DigitsB~147 spike falls into.
BOOST_AUTO_TEST_CASE(estimate_reflects_genuine_near_singularity)
{
	using complex_mp = bertini::complex_mp;
	using real_mp   = bertini::real_mp;

	for (int k = 20; k <= 140; k += 40)
	{
		bertini::DefaultPrecision(static_cast<unsigned int>(k + 50)); // enough digits to represent eps = 1e-k

		real_mp eps = pow(real_mp(10), -k);
		// J = [[1, 1], [1, 1+eps]] : det = eps, so ||J^{-1}||_2 ~ 2/eps ~ 1e+k (analytic).
		Mat<complex_mp> J(2, 2);
		J(0,0) = complex_mp(1); J(0,1) = complex_mp(1);
		J(1,0) = complex_mp(1); J(1,1) = complex_mp(1) + complex_mp(eps);

		auto lu = J.partialPivLu();
		real_mp best(0);
		for (int i = 0; i < 12; ++i)
		{
			Vec<complex_mp> r = bertini::RandomOfUnits<complex_mp>(2);
			real_mp nrm = lu.solve(r).norm();
			if (nrm > best) best = nrm;
		}

		double log_est = static_cast<double>(log10(best));
		double log_tru = static_cast<double>(log10(real_mp(2) / eps)); // ~ k

		// the estimate is the right order of magnitude (within ~2 decades) of the true value --
		// so a DigitsB spike of ~147 means a genuinely ~1e145 ||J^{-1}||, NOT an estimation artifact.
		BOOST_CHECK_SMALL(log_est - log_tru, 2.0);   // est <= tru * 1e2
		BOOST_CHECK_GT(log_est, log_tru - 2.0);      // est >= tru / 1e2
	}
}


BOOST_AUTO_TEST_SUITE_END()
