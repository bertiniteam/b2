//This file is part of Bertini 2.
//
//total_degree_linear_product.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//total_degree_linear_product.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with total_degree_linear_product.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file total_degree_linear_product.hpp

\brief Defines the TotalDegreeLinearProduct start system type.
*/

#pragma once

#include "bertini2/system/start_base.hpp"
#include "bertini2/system/start/utility.hpp"

namespace bertini
{
	namespace start_system{


		/**
		\brief Total degree start system for 1-homogeneous polynomial systems, built from random linear products.

		The standard, well-conditioned total-degree start system in Numerical Algebraic Geometry.

		For a square target system with a single affine variable group of \f$n\f$ variables, the
		start function for a degree-\f$d_i\f$ target function is a product of \f$d_i\f$ random affine
		linear forms, \f$\prod_{j=1}^{d_i} (a_{ij}\cdot x + b_{ij})\f$.  Its start points are generic
		intersections (in general position), found by linear algebra -- NOT a roots-of-unity lattice.
		The number of start points is the Bezout bound \f$\prod_i d_i\f$.  For the (cheaper, structured)
		roots-of-unity start, see TotalDegreeBinomial.

		This is the single-affine-variable-group specialization of MHomogeneous: like MHom, the start
		system evaluates through a products-of-linears block (the SLP compiler cannot compile
		linear-product node trees), and each start point is the solution of an \f$n\times n\f$ linear
		system formed from one chosen linear factor per function.

		Note that the corresponding target system MUST be square -- have the same number of functions
		and variables.  The start system cannot be constructed otherwise; it throws if not square.

		The start points are accessed by index (unsigned long long), instead of being generated all at once.
		*/
		class TotalDegreeLinearProduct : public StartSystem
		{
		public:
			TotalDegreeLinearProduct() = default;
			virtual ~TotalDegreeLinearProduct() = default;

			/**
			 Constructor for making a total degree start system from a polynomial system

			 \throws std::runtime_error, if the input target system is not square, is not polynomial, has a path variable already, has more than one variable group, or has any homogeneous variable groups.
			*/
			TotalDegreeLinearProduct(System const& s);


			/**
			Get the number of start points for this total degree start system.  This is the Bezout bound for the target system.  Provided here for your convenience.
			*/
			unsigned long long NumStartPoints() const override;

			/// \brief Multiply this start system in place by an expression node.
			TotalDegreeLinearProduct& operator*=(Nd const& n);

			TotalDegreeLinearProduct& operator+=(System const& sys) = delete;

			/// \brief Check that the target system is suitable for a total-degree start system.
			void SanityChecks(System const& s);

		private:

			/**
			Copy the degrees from another system into this one
			*/
			void CopyDegrees(System const& s);

			/**
			Populate the random linear-factor coefficients of this start system.  For each function i
			(degree d_i) this fills linear_coeffs_[i], a (d_i) x (n+1) matrix whose rows are the random
			affine linear factors (the trailing column being each factor's constant term).
			*/
			void SeedLinearCoeffs(System const& s);

			/**
			Build the products-of-linears evaluation block from linear_coeffs_ (so the start system
			evaluates via the block, as MHomogeneous does).  Assumes the variable structure is set up
			(homogenized/patched as appropriate) and the coefficients are seeded.
			*/
			void BuildBlock(System const& s);


			/**
			Get the ith start point, in double precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<complex_dbl> GenerateStartPoint(complex_dbl,unsigned long long index) const override;

			/**
			Get the ith start point, in current default precision.

			Called by the base StartSystem's StartPoint(index) method.
			*/
			Vec<complex_mp> GenerateStartPoint(complex_mp,unsigned long long index) const override;

			/**
			 A local version of GenerateStartPoint that can be templated.
			*/
			template<typename T>
			void GenerateStartPointT(Vec<T>& start_point, unsigned long long index) const;

			std::vector<unsigned long long> degrees_; ///< stores the degrees of the functions.
			/// Random linear-factor coefficients, one matrix per function.  Entry i is a
			/// (degrees_[i]) x (NumNaturalVariables()+1) matrix: each row is an affine linear factor
			/// over the variables, the trailing column being the factor's constant term.  Generated at
			/// MaxPrecisionAllowed so the block's master is precision-faithful.  This is the
			/// single-group analogue of MHomogeneous::linear_coeffs_.
			std::vector<Mat<complex_mp>> linear_coeffs_;


			friend class boost::serialization::access;

			template <typename Archive>
			void serialize(Archive& ar, const unsigned /*version*/) {
				ar & boost::serialization::base_object<StartSystem>(*this);
				ar & degrees_;
				// serialize the coefficients too (unlike MHomogeneous, which persists only degrees_):
				// a round-tripped TotalDegreeLinearProduct can then regenerate its start points, not merely evaluate.
				ar & linear_coeffs_;
			}

		};
	}
}
