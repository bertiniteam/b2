//This file is part of Bertini 2.
//
//slice.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//slice.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with slice.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file bertini2/system/slice.hpp

\brief Provides the bertini::Slice class -- a linear slice of affine or projective space.

A Slice is the linear part of a witness set: a stack of linear forms whose common zero set
cuts a positive-dimensional component down to isolated witness points.  It is backed by a
bertini::blocks::LinearFormsBlock (the augmented coefficient matrix M, one row per linear
form, the last column carrying each form's constant term), so a Slice evaluates, carries
precision, and drops into a System exactly like any other evaluation block.  Because that
augmented (num_vars+1)-column layout is shared with ProductsOfLinearsBlock, a Slice's rows
are also ready-made factors for a product-of-linears block -- the composition regeneration
needs.

\see bertini::blocks::LinearFormsBlock
\see bertini::blocks::ProductsOfLinearsBlock
*/


#ifndef BERTINI_SLICE_HPP
#define BERTINI_SLICE_HPP

#include <vector>

#include "bertini2/function_tree.hpp"
#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/system/blocks/linear_forms_block.hpp"


namespace bertini {

	// Slice::AddTo hands the slice's linear-forms block to a System; we only need System's
	// name here (the definition lives in slice.cpp, which includes system.hpp).
	class System;

	/**
	\brief A linear slice of affine or projective space: a stack of linear forms, M [x ; 1].

	The slice is held as an augmented coefficient matrix (one row per linear form, the
	trailing column being that form's constant term) inside a LinearFormsBlock.  A
	homogeneous slice simply has a zero constant column.  Slices compose: Head / Tail / Rows
	return a new Slice over the same variables built from a subset of the linear forms.
	*/
	class Slice
	{
		blocks::LinearFormsBlock block_;  ///< the linear forms M [x ; 1] -- eval / precision / coefficients
		VariableGroup sliced_vars_;       ///< the variables this slice is a function of
		bool is_homogeneous_ = false;     ///< whether the forms were authored without constant terms

	public:

		/// An empty slice (zero forms, zero variables).
		Slice() = default;

		/**
		\brief Build a slice directly from an augmented coefficient matrix.

		\param v The variables the slice is a function of.
		\param augmented_coefficients One row per linear form, (number-of-forms) x (v.size()+1);
		       the trailing column is each form's constant term (zero, for a homogeneous slice).
		\param homogeneous Whether the slice was authored without constant terms.
		*/
		static Slice FromCoefficients(VariableGroup const& v, Mat<complex_mp> const& augmented_coefficients, bool homogeneous = false)
		{
			assert(static_cast<size_t>(augmented_coefficients.cols()) == v.size() + 1 &&
			       "a slice coefficient matrix must have (num_variables + 1) columns");
			Slice s;
			s.sliced_vars_ = v;
			s.is_homogeneous_ = homogeneous;
			s.block_ = blocks::LinearFormsBlock(v.size(), augmented_coefficients);
			return s;
		}

		/**
		\brief Produce a random real slice on a variable group, slicing a given number of dimensions.
		*/
		static Slice RandomReal(VariableGroup const& v, unsigned dim, bool homogeneous = false, bool orthogonal = true)
		{
			typedef void (*funtype) (complex_mp&, unsigned); // the type for number generation
			// bounded-modulus draw (away from 0 and infinity), matching patches and the start systems;
			// kept REAL so a real slice stays real.  (A deeper pass on slice generation -- the constant
			// column and the orthogonal=false path -- is still TODO.)
			funtype gen = bertini::multiprecision::RandomRealBoundedModulusAssign;
			return Make(v, dim, homogeneous, orthogonal, gen);
		}

		/**
		\brief Generate a random complex slice.
		*/
		static Slice RandomComplex(VariableGroup const& v, unsigned dim, bool homogeneous = false, bool orthogonal = true)
		{
			typedef void (*funtype) (complex_mp&, unsigned); // the type for number generation
			// bounded-modulus draw (away from 0 and infinity), matching patches and the start systems.
			funtype gen = bertini::multiprecision::RandomComplexBoundedModulusAssign;
			return Make(v, dim, homogeneous, orthogonal, gen);
		}

		/**
		\brief Factory for generating slices.  Generates the variable-coefficient block (optionally
		orthonormalized by a QR factorization) and the constant column, then assembles the augmented
		matrix the LinearFormsBlock holds.
		*/
		static Slice Make(VariableGroup const& v, unsigned dim, bool homogeneous, bool orthogonal, std::function<void(complex_mp&, unsigned)> gen)
		{
			const unsigned num_vars = static_cast<unsigned>(v.size());

			Mat<complex_mp> coeffs(dim, num_vars); // the variable coefficients (one row per form)

			if (orthogonal)
			{
				using std::min;
				using std::max;

				auto mindim = min(dim, num_vars);
				auto maxdim = max(dim, num_vars);

				bool need_transpose = dim < num_vars;

				coeffs.resize(maxdim, mindim);

				for (unsigned ii(0); ii < maxdim; ++ii)
					for (unsigned jj(0); jj < mindim; ++jj)
						gen(coeffs(ii, jj), MaxPrecisionAllowed());

				auto prev_precision = DefaultPrecision();
				DefaultPrecision(MaxPrecisionAllowed());

				auto QR_factorization = Eigen::HouseholderQR<Mat<complex_mp> >(coeffs);
				coeffs = QR_factorization.householderQ() * Mat<complex_mp>::Identity(maxdim, mindim);

				if (need_transpose)
					coeffs.transposeInPlace();

				DefaultPrecision(prev_precision);
			}
			else
			{
				for (unsigned ii(0); ii < dim; ++ii)
					for (unsigned jj(0); jj < num_vars; ++jj)
						gen(coeffs(ii, jj), MaxPrecisionAllowed());
			}

			assert(static_cast<unsigned>(coeffs.rows()) == dim);
			assert(static_cast<unsigned>(coeffs.cols()) == num_vars);

			// Assemble the augmented matrix: [ coeffs | constants ].  A homogeneous slice's constant
			// column is zero; otherwise it is freshly generated.
			Mat<complex_mp> augmented(dim, num_vars + 1);
			augmented.leftCols(num_vars) = coeffs;
			if (homogeneous)
				augmented.col(num_vars).setZero();
			else
				for (unsigned ii(0); ii < dim; ++ii)
					gen(augmented(ii, num_vars), MaxPrecisionAllowed());

			return FromCoefficients(v, augmented, homogeneous);
		}


		/**
		\brief Evaluate the slice's linear-form values, in-place.
		*/
		template<typename NumT>
		void Eval(Vec<NumT> & result, Vec<NumT> const& x) const
		{
			result.resize(Dimension());
			NumT path_value(0); // ignored: linear forms are autonomous
			block_.EvalInPlace<NumT>(result, x, path_value);
		}

		/**
		\brief Evaluate the slice's linear-form values.
		*/
		template<typename NumT>
		Vec<NumT> Eval(Vec<NumT> const& x) const
		{
			Vec<NumT> result(Dimension());
			Eval(result, x);
			return result;
		}

		/**
		\brief The slice's Jacobian (its constant variable-coefficient matrix), in-place.
		*/
		template<typename NumT>
		void Jacobian(Mat<NumT> & result, Vec<NumT> const& x) const
		{
			result.resize(Dimension(), NumVariables());
			NumT path_value(0);
			block_.JacobianInPlace<NumT>(result, x, path_value);
		}

		/**
		\brief The slice's Jacobian (its constant variable-coefficient matrix).
		*/
		template<typename NumT>
		Mat<NumT> Jacobian(Vec<NumT> const& x) const
		{
			Mat<NumT> result(Dimension(), NumVariables());
			Jacobian(result, x);
			return result;
		}


		/**
		\brief The augmented coefficient matrix: one row per linear form, (Dimension) x (NumVariables+1),
		the trailing column carrying each form's constant term.

		These rows are also factor rows for a ProductsOfLinearsBlock, so a slice composes directly into
		the product-of-linears form regeneration uses.
		*/
		Mat<complex_mp> const& Coefficients() const
		{
			return block_.Coefficients();
		}

		/// The underlying linear-forms block (eval / Jacobian / precision engine).
		blocks::LinearFormsBlock const& AsLinearFormsBlock() const
		{
			return block_;
		}

		/// Add this slice's linear forms to a System as a LinearFormsBlock.  (Defined in slice.cpp.)
		void AddTo(System & s) const;

		/// A standalone System whose functions are exactly this slice's linear forms (over the slice's
		/// variable group).  Lets a slice be carried around and evaluated / tracked on its own.
		/// (Defined in slice.cpp.)
		System AsSystem() const;

		/**
		\brief A new slice stacking this slice's linear forms on top of \p other's.

		Both slices must be on the same number of variables.  The result is homogeneous only if both
		operands are.  This is how you build a higher-codimension slice from pieces (and the Python
		`+` operator).
		*/
		Slice Concatenate(Slice const& other) const
		{
			if (NumVariables() != other.NumVariables())
				throw std::runtime_error("Slice::Concatenate requires both slices to be on the same number of variables");

			Mat<complex_mp> const& A = Coefficients();
			Mat<complex_mp> const& B = other.Coefficients();
			Mat<complex_mp> stacked(A.rows() + B.rows(), A.cols());
			stacked.topRows(A.rows()) = A;
			stacked.bottomRows(B.rows()) = B;

			return FromCoefficients(sliced_vars_, stacked, is_homogeneous_ && other.is_homogeneous_);
		}


		/**
		\brief A new slice over the same variables built from the first \p m linear forms.
		*/
		Slice Head(unsigned m) const
		{
			if (m > Dimension())
				throw std::runtime_error("Slice::Head asked for more forms than the slice has");
			return FromCoefficients(sliced_vars_, Coefficients().topRows(m), is_homogeneous_);
		}

		/**
		\brief A new slice over the same variables built from the last \p m linear forms.
		*/
		Slice Tail(unsigned m) const
		{
			if (m > Dimension())
				throw std::runtime_error("Slice::Tail asked for more forms than the slice has");
			return FromCoefficients(sliced_vars_, Coefficients().bottomRows(m), is_homogeneous_);
		}

		/**
		\brief A new slice over the same variables built from the chosen linear forms.
		*/
		Slice Rows(std::vector<unsigned> const& indices) const
		{
			Mat<complex_mp> const& C = Coefficients();
			Mat<complex_mp> sub(static_cast<Eigen::Index>(indices.size()), C.cols());
			for (size_t ii = 0; ii < indices.size(); ++ii)
			{
				if (indices[ii] >= Dimension())
					throw std::runtime_error("Slice::Rows asked for a form index outside the slice");
				sub.row(static_cast<Eigen::Index>(ii)) = C.row(indices[ii]);
			}
			return FromCoefficients(sliced_vars_, sub, is_homogeneous_);
		}


		/**
		\brief The dimension of the slice -- the number of linear forms.
		*/
		unsigned Dimension() const
		{
			return static_cast<unsigned>(block_.NumFunctions());
		}

		/**
		\brief The number of variables sliced.
		*/
		unsigned NumVariables() const
		{
			return static_cast<unsigned>(sliced_vars_.size());
		}

		/**
		\brief The variables the slice is a function of.
		*/
		VariableGroup const& Variables() const
		{
			return sliced_vars_;
		}

		/**
		\brief Whether the slice was authored without constant terms (passes through the origin).
		*/
		bool IsHomogeneous() const
		{
			return is_homogeneous_;
		}


		/**
		\brief Get the current working precision of the slice, in digits.
		*/
		unsigned Precision() const
		{
			return block_.Precision();
		}

		/**
		\brief Set the working precision of the slice, in digits.
		*/
		void Precision(unsigned new_precision) const
		{
			block_.Precision(new_precision);
		}

	private:

		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & block_;
			ar & sliced_vars_;
			ar & is_homogeneous_;
		}

		friend std::ostream& operator<<(std::ostream&, Slice const&);
	};

	/**
	\brief Provides output streaming for Slice
	*/
	std::ostream& operator<<(std::ostream& out, Slice const& s);
} // re: namespace bertini

#endif
