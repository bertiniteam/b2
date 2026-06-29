//This file is part of Bertini 2.
//
//total_degree.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//total_degree.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with total_degree.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

#include "bertini2/system/start/total_degree.hpp"

#include "bertini2/system/blocks/block.hpp"
#include "bertini2/random.hpp"

#include <map>


BOOST_CLASS_EXPORT(bertini::start_system::TotalDegree);


namespace bertini {
	using namespace bertini::node;

	namespace start_system {

		// constructor for TotalDegree start system, from any other *suitable* system.
		// This is the single-affine-variable-group specialization of MHomogeneous: one variable
		// group, no homogeneous groups, no partitions.  See mhom.cpp for the multi-group version.
		TotalDegree::TotalDegree(System const& s)
		{
			SanityChecks(s);
			CopyDegrees(s);
			CopyVariableStructure(s);
			SeedLinearCoeffs(s);

			if (s.IsHomogeneous())
				Homogenize();

			if (s.IsPatched())
				CopyPatches(s);

			BuildBlock(s);
		}// total degree constructor


		TotalDegree& TotalDegree::operator*=(Nd const& n)
		{
			System::operator*=(n);
			return *this;
		}


		unsigned long long TotalDegree::NumStartPoints() const
		{
			unsigned long long num_start_points = 1;
			for (const auto& iter : degrees_)
				num_start_points*=iter;
			return num_start_points;
		}


		void TotalDegree::SanityChecks(System const& s)
		{
			if (s.NumHomVariableGroups() > 0)
				throw std::runtime_error("a homogeneous variable group is present.  currently unallowed");

			if (s.NumTotalFunctions() != s.NumVariables())
				throw std::runtime_error("attempting to construct total degree start system from non-square target system");

			if (s.HavePathVariable())
				throw std::runtime_error("attempting to construct total degree start system, but target system has path varible declared already");

			if (s.NumVariableGroups() != 1)
				throw std::runtime_error("more than one affine variable group.  currently unallowed");

			if (!s.IsPolynomial())
				throw std::runtime_error("attempting to construct total degree start system from non-polynomial target system");
		}


		void TotalDegree::CopyDegrees(System const& s)
		{
			auto deg = s.Degrees();
			for (const auto& d : deg)
				degrees_.push_back(static_cast<size_t>(d));
		}


		// Generate the random linear-factor coefficients directly (no node::LinearProduct), exactly
		// as mhom.cpp does for one (function, group) pair.  For each function i (degree d_i) fill a
		// d_i x (n+1) matrix: each row is an affine linear factor over the n natural variables, the
		// trailing column being the factor's constant.  Generated at MaxPrecisionAllowed so the
		// block's master is precision-faithful.
		void TotalDegree::SeedLinearCoeffs(System const& s)
		{
			auto const saved_prec = DefaultPrecision();
			DefaultPrecision(MaxPrecisionAllowed());

			const Eigen::Index n = static_cast<Eigen::Index>(s.NumNaturalVariables());
			linear_coeffs_.clear();
			linear_coeffs_.reserve(degrees_.size());
			for (size_t ii = 0; ii < degrees_.size(); ++ii)
			{
				const Eigen::Index d = static_cast<Eigen::Index>(degrees_[ii]);
				Mat<complex_mp> C(d, n + 1);
				for (Eigen::Index f = 0; f < d; ++f)
					for (Eigen::Index k = 0; k < n + 1; ++k)   // includes the trailing constant column
						C(f, k) = complex_mp(real_mp(RandomRat()), real_mp(RandomRat()));
				linear_coeffs_.push_back(std::move(C));
			}

			DefaultPrecision(saved_prec);
		}


		// Build the products-of-linears evaluation block from linear_coeffs_, so the start system
		// evaluates via the block (the SLP compiler cannot compile linear-product node trees).  Per
		// function we assemble an augmented coefficient matrix (one row per factor, length
		// NumVariables()+1): each factor's coefficients are placed by VARIABLE IDENTITY into the full
		// variable ordering, and the factor's constant multiplies the group's homogenizing variable
		// when the system is homogeneous (so it goes in that variable's column), else it is the
		// augmented constant in the trailing column.  This is mhom.cpp's block assembly with a single
		// affine group (no projective groups).  Built at MaxPrecisionAllowed.
		void TotalDegree::BuildBlock(System const& /*s*/)
		{
			auto const saved_prec = DefaultPrecision();
			DefaultPrecision(MaxPrecisionAllowed());
			this->precision(MaxPrecisionAllowed());

			const VariableGroup& vars = this->Variables();
			std::map<node::Node const*, Eigen::Index> col_of;
			for (Eigen::Index c = 0; c < static_cast<Eigen::Index>(vars.size()); ++c)
				col_of[vars[static_cast<size_t>(c)].get()] = c;
			const Eigen::Index n = static_cast<Eigen::Index>(this->NumVariables());

			// the single affine group's variables, and (if homogenized) its homogenizing variable.
			const VariableGroup& gvars = this->AffineVariableGroup(0);
			std::shared_ptr<node::Variable> hom_var;
			if (!this->HomogenizingVariables().empty())
				hom_var = this->HomogenizingVariables()[0];

			std::vector<Mat<complex_mp>> per_function;
			per_function.reserve(linear_coeffs_.size());

			for (size_t ii = 0; ii < linear_coeffs_.size(); ++ii)
			{
				const Mat<complex_mp>& C = linear_coeffs_[ii];
				const Eigen::Index d = C.rows();
				Mat<complex_mp> M(d, n + 1);
				for (Eigen::Index f = 0; f < d; ++f)
				{
					Vec<complex_mp> row = Vec<complex_mp>::Zero(n + 1);
					for (size_t k = 0; k < gvars.size(); ++k)
						row(col_of.at(gvars[k].get())) = C(f, static_cast<Eigen::Index>(k));
					const complex_mp& constant = C(f, static_cast<Eigen::Index>(gvars.size()));
					if (hom_var && col_of.count(hom_var.get()))
						row(col_of.at(hom_var.get())) = constant;
					else
						row(n) = constant;   // augmented constant (affine, un-homogenized)
					M.row(f) = row.transpose();
				}
				per_function.push_back(std::move(M));
			}

			this->AddBlock(blocks::ProductsOfLinearsBlock(static_cast<size_t>(n), std::move(per_function)));

			DefaultPrecision(saved_prec);
			this->precision(saved_prec);
		}


		template<typename T>
		void TotalDegree::GenerateStartPointT(Vec<T>& start_point, unsigned long long index) const
		{
			// Decompose the flat index into one chosen linear factor per function.  (Single group =>
			// no partition search, unlike MHomogeneous.)
			std::vector<size_t> dim_vector(degrees_.size());
			for (size_t ii = 0; ii < degrees_.size(); ++ii)
				dim_vector[ii] = static_cast<size_t>(degrees_[ii]);
			std::vector<size_t> subscript = IndexToSubscript<size_t>(index, dim_vector);

			// Build the n x n linear system whose solution is this start point, in natural
			// coordinates: function i contributes its chosen affine linear factor = 0, i.e. row i is
			// the factor's variable coefficients and the right-hand side is minus its constant.
			const Eigen::Index n = static_cast<Eigen::Index>(NumNaturalVariables());
			Mat<T> A = Mat<T>::Zero(n, n);
			Vec<T> b = Vec<T>::Zero(n);

			// coefficients come from linear_coeffs_ (mpfr master); cast each to the working type T
			// (a no-op widen for mpfr, a narrowing for complex_dbl).
			auto as_T = [](complex_mp const& z) -> T {
				if constexpr (std::is_same<T, complex_dbl>::value) return complex_dbl(z);
				else return z;
			};

			for (Eigen::Index ii = 0; ii < n; ++ii)
			{
				const Mat<complex_mp>& C = linear_coeffs_[static_cast<size_t>(ii)];
				const Eigen::Index f = static_cast<Eigen::Index>(subscript[static_cast<size_t>(ii)]);
				for (Eigen::Index k = 0; k < n; ++k)
					A(ii, k) = as_T(C(f, k));
				b(ii) = -as_T(C(f, n));   // move the constant term to the right-hand side
			}

			Vec<T> affine_solution = A.partialPivLu().solve(b);

			// The linear solve gives the start point in affine (dehomogenized) coordinates; lift it
			// onto the homogenized + patched coordinate system the homotopy is tracked in
			// (HomogenizePoint inserts the homogenizing coordinate and rescales onto the patch, and is
			// a no-op if not homogenized).
			start_point = this->HomogenizePoint(affine_solution);

			// Return the point at the current working precision.  The coefficients and the patch are
			// stored at MaxPrecisionAllowed, which would otherwise leak into the start point and
			// mismatch the tracker's working precision (the adaptive tracker begins at the ambient
			// precision).
			if constexpr (!std::is_same<T, complex_dbl>::value)
				for (Eigen::Index i = 0; i < start_point.size(); ++i)
					start_point(i).precision(DefaultPrecision());
		}


		Vec<complex_dbl> TotalDegree::GenerateStartPoint(complex_dbl,unsigned long long index) const
		{
			Vec<complex_dbl> start_point(NumVariables());
			GenerateStartPointT(start_point, index);
			return start_point;
		}


		Vec<complex_mp> TotalDegree::GenerateStartPoint(complex_mp,unsigned long long index) const
		{
			Vec<complex_mp> start_point(NumVariables());
			GenerateStartPointT(start_point, index);
			return start_point;
		}

		inline
		TotalDegree operator*(TotalDegree td, std::shared_ptr<node::Node> const& n)
		{
			td *= n;
			return td;
		}
	} // namespace start_system
} //namespace bertini
