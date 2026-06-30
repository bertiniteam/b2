//This file is part of Bertini 2.
//
//eigen_extensions.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//eigen_extensions.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with eigen_extensions.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, University of Wisconsin Eau Claire

/**
\file eigen_extensions.hpp

\brief Bertini extensions of the Eigen linear algebra library.
*/

#ifndef BERTINI_EIGEN_EXTENSIONS_HPP
#define BERTINI_EIGEN_EXTENSIONS_HPP


#include "bertini2/num_traits.hpp"

#include <boost/serialization/complex.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/split_member.hpp>

/// \brief Tell Eigen to mix a Boost.serialization addon into every dense matrix/array.
#define EIGEN_DENSEBASE_PLUGIN "bertini2/eigen_serialization_addon.hpp"

#include <Eigen/Core>

namespace {
using mpfr_real = bertini::real_mp;
using complex_mp = bertini::complex_mp;
}

namespace Eigen {

/// \cond EIGEN_GLUE
	template<> struct NumTraits<mpfr_real> : GenericNumTraits<mpfr_real> // permits to get the epsilon, dummy_precision, lowest, highest functions
	{

		typedef mpfr_real Real;
		typedef mpfr_real NonInteger;
		typedef mpfr_real Nested;
		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned = 1,
			RequireInitialization = 1, // yes, require initialization, otherwise get crashes
			ReadCost = 20,
			AddCost = 30,
			MulCost = 40
		};


		inline static Real highest() {

			return (mpfr_real(1) - epsilon()) * pow(mpfr_real(2),mpfr_get_emax()-1);//);//DefaultPrecision());
		}

		inline static Real lowest() {
			return -highest();
		}

		// ThreadPrecision (thread-local), not DefaultPrecision (global): these
		// tolerances are consumed inside tracking (LU solves, norms, convergence
		// checks), which may run on std::thread workers whose precision is set
		// thread-locally.  On the main thread the two agree.
		// epsilon() and dummy_precision() are pow(10, -precision): expensive (transcendental,
		// heap-allocating) yet a pure function of the thread precision.  They are called per
		// LU-diagonal element, per Newton iteration, per predictor stage, per step, so recomputing
		// them dominated the multiprecision path.  Memoize per thread, keyed on the current
		// precision -- the cached value is bit-identical to a fresh pow at that precision.
		inline static Real dummy_precision()
		{
			using bertini::ThreadPrecision;
			thread_local unsigned cached_prec = 0;
			thread_local Real cached_val;
			const unsigned p = ThreadPrecision();
			if (p != cached_prec)
			{
				cached_val = pow(mpfr_real(10), -int(p-3));
				cached_prec = p;
			}
			return cached_val;
		}

		inline static Real epsilon()
		{
			using bertini::ThreadPrecision;
			thread_local unsigned cached_prec = 0;
			thread_local Real cached_val;
			const unsigned p = ThreadPrecision();
			if (p != cached_prec)
			{
				cached_val = pow(mpfr_real(10), -int(p));
				cached_prec = p;
			}
			return cached_val;
		}

		static inline int digits10()
		{
			return static_cast<int>(bertini::ThreadPrecision());
			// return internal::default_digits10_impl<T>::run();
		}
		//http://www.manpagez.com/info/mpfr/mpfr-2.3.2/mpfr_31.php
	};


	template<typename Expr1,typename Expr2, typename Expr3>
	struct NumTraits<
		boost::multiprecision::detail::expression<Expr1,
      	Expr2, Expr3, void, void>
      		> : NumTraits<mpfr_real> // permits to get the epsilon, dummy_precision, lowest, highest functions
	{

		typedef mpfr_real Real;
		typedef mpfr_real NonInteger;
		typedef mpfr_real Nested;

		//http://www.manpagez.com/info/mpfr/mpfr-2.3.2/mpfr_31.php
	};



	/**
	 \brief This templated struct permits us to use the complex_mp type in Eigen matrices.

	 Provides methods to get the epsilon, dummy_precision, lowest, highest functions, largely by inheritance from the NumTraits<mpfr_real> contained in mpfr_extensions.
	 */
	template<> struct NumTraits<complex_mp> : NumTraits<mpfr_real>
	{
		typedef mpfr_real Real;
		typedef mpfr_real NonInteger;
		typedef complex_mp Nested;// Nested;
		enum {
			IsComplex = 1,
			IsInteger = 0,
			IsSigned = 1,
			RequireInitialization = 1, // yes, require initialization, otherwise get crashes
			ReadCost = 2 * NumTraits<Real>::ReadCost,
			AddCost = 2 * NumTraits<Real>::AddCost,
			MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
		};

	};

	namespace internal {
		template<>
		struct abs2_impl<complex_mp>
		{
			static inline mpfr_real run(const complex_mp& x)
			{
				return real(x)*real(x) + imag(x)*imag(x);
			}
		};


		template<> inline complex_mp random<complex_mp>()
		{
			return bertini::multiprecision::rand();
		}

		template<> inline complex_mp random<complex_mp>(const complex_mp& a, const complex_mp& b)
		{
			return a + (b-a) * random<complex_mp>();
		}

		template<>
		struct conj_helper<complex_mp, complex_mp, false, true>
		{
			typedef complex_mp Scalar;
			EIGEN_STRONG_INLINE Scalar pmadd(const Scalar& x, const Scalar& y, const Scalar& c) const
			{ return c + pmul(x,y); }

			EIGEN_STRONG_INLINE Scalar pmul(const Scalar& x, const Scalar& y) const
			{ return Scalar(real(x)*real(y) + imag(x)*imag(y), imag(x)*real(y) - real(x)*imag(y)); }
		};

		template<>
		struct conj_helper<complex_mp, complex_mp, true, false>
		{
			typedef complex_mp Scalar;
			EIGEN_STRONG_INLINE Scalar pmadd(const Scalar& x, const Scalar& y, const Scalar& c) const
			{ return c + pmul(x,y); }

			EIGEN_STRONG_INLINE Scalar pmul(const Scalar& x, const Scalar& y) const
			{ return Scalar(real(x)*real(y) + imag(x)*imag(y), real(x)*imag(y) - imag(x)*real(y)); }
		};

		// see https://forum.kde.org/viewtopic.php?f=74&t=111176

		//int
		template<>
		struct scalar_product_traits<int,complex_mp>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};

		template<>
		struct scalar_product_traits<complex_mp, int>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};

		//long
		template<>
		struct scalar_product_traits<long,complex_mp>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};

		template<>
		struct scalar_product_traits<complex_mp, long>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};

		//long long
		template<>
		struct scalar_product_traits<long long,complex_mp>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};

		template<>
		struct scalar_product_traits<complex_mp, long long>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};


		//mpfr_real
		template<>
		struct scalar_product_traits<mpfr_real,complex_mp>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};

		template<>
		struct scalar_product_traits<complex_mp, mpfr_real>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};

		//mpz_int
		template<>
		struct scalar_product_traits<bertini::mpz_int,complex_mp>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};

		template<>
		struct scalar_product_traits<complex_mp, bertini::mpz_int>
		{
	    	enum { Defined = 1 };
	    	typedef complex_mp ReturnType;
		};

	} // re: namespace internal
/// \endcond
} // re: namespace Eigen


#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <Eigen/Dense>

namespace bertini {

	/// \brief A dynamically-sized column vector of NumT.
	template<typename NumT> using Vec = Eigen::Matrix<NumT, Eigen::Dynamic, 1>;
	/// \brief A dynamically-sized matrix of NumT.
	template<typename NumT> using Mat = Eigen::Matrix<NumT, Eigen::Dynamic, Eigen::Dynamic>;


	/**
	\brief Checks whether the number of rows and columns are either 0
	*/
	template<typename Derived>
	inline
	bool IsEmpty(Eigen::MatrixBase<Derived> const & v)
	{
		return (v.rows()==0) || (v.cols()==0);
	}


	/**
	\brief Get the precision of an Eigen object. If the object is empty, it's the precision of a default-constructed Scalar.  If it actually has content, then it's the precision of the first element.

	If you require that the object be of uniform precision when you check, use PrecisionRequireUniform
	*/
	template<typename Derived>
	inline
	unsigned Precision(Eigen::MatrixBase<Derived> const & v)
	{
		if (IsEmpty(v))
			throw std::runtime_error("getting precision of empty object");

		return Precision(v(0,0));
	}


	/**
	\brief Query the precision of an Eigen object, requiring it to be non-empty, and of uniform precision

	\throw runtime_error if not uniform precision
	Dies in an Eigen assert if empty (assuming you didn't disable these)
	*/
	template<typename Derived>
	unsigned PrecisionRequireUniform(Eigen::MatrixBase<Derived> const & v)
	{
		const auto a = Precision(v(0,0));
		for (int ii=0; ii<v.cols(); ++ii)
			for (int jj=0; jj<v.rows(); ++jj)
			{
				if (a!=Precision(v(jj,ii)))
				{
					std::stringstream ss;
					ss << "non-uniform precision in object! (" << a << "!=" << Precision(v(jj,ii)) <<   ") at position " << jj << "," << ii;
					throw std::runtime_error(ss.str());
				}
			}
		return a;
	}


	/**
	\brief Set the precision of an Eigen object.

	\param v the matrix to change
	\param prec the new precision
	*/
	template<typename Derived>
	void Precision(Eigen::MatrixBase<Derived> & v, unsigned prec)
	{
		using bertini::Precision;
		for (int ii=0; ii<v.cols(); ++ii)
			for (int jj=0; jj<v.rows(); ++jj)
				Precision(v(jj,ii),prec);
	}



	/**
	 Test a numbers being very small.

	 Compares the number against machine epsilon (or software epsilon if a multiple precision type), times 100.

	 \note Machine epsilon for doubles is about 1e-16, for mpfr_real, it's 10^-current precision.

	 \param testme The number to test

	 \return true, if the number is very small.  False otherwise.
	*/
	template<typename T>
	inline
	bool IsSmallValue(T const& testme)
	{
		// |testme| <= eps*100  <=>  abs2(testme) <= (eps*100)^2.  Comparing squared magnitudes uses
		// abs2 (re^2+im^2 for complex; see abs2_impl<complex_mp>) and avoids the allocating, sqrt-based
		// std::abs/hypot on the multiprecision complex type.  eps*100 is non-negative, so the squared
		// comparison is exactly equivalent.
		using Real = typename Eigen::NumTraits<T>::Real;
		const Real thresh = Eigen::NumTraits<T>::epsilon() * 100;
		return Eigen::numext::abs2(testme) <= thresh * thresh;
	}

	/**
	\brief Check whether two values are very close to each other.

	\tparam T The numeric type.
	\param a The first value.
	\param b The second value.
	\param e The tolerance for being close.

	See http://www.boost.org/doc/libs/1_34_0/libs/test/doc/components/test_tools/floating_point_comparison.html, for example.
	*/
	template<typename T>
	inline
	bool IsSymmRelDiffSmall(T const& a, T const& b, typename Eigen::NumTraits<T>::Real const& e)
	{
		using std::abs;
		if (a==b)
			return true;

		typename Eigen::NumTraits<T>::Real c = abs(a-b);
		return (c/abs(a) <= e) || (c/abs(b) <= e) ;
	}

	/**
	 Test two numbers for having large ratio.

	 The basis for comparison is Eigen's fuzzy precision, Eigen::NumTraits<T>::dummy_precision();

	 \note For doubles, the threshold is 1e-12, for mpfr_real is 1e3*current precision.

	 \param numerator The numerator for the ratio.
	 \param denomenator The denomenator for the ratio.

	 \return true, if the ratio is very large, false otherwise
	*/
	template<typename T>
	inline
	bool IsLargeChange(T const& numerator, T const& denomenator)
	{
		static_assert(!Eigen::NumTraits<T>::IsInteger, "IsLargeChange cannot be used safely on non-integral types");
		// |num/den| >= 1/dummy  <=>  |num|*dummy >= |den|  <=>  abs2(num)*dummy^2 >= abs2(den).
		// The squared form drops BOTH the complex division (num/den) and the sqrt-based abs/hypot --
		// each of which allocates multiprecision temporaries -- and is exact since all quantities are
		// non-negative.  (In its hot caller, LUPartialPivotDecompositionSuccessful, the denominator has
		// already passed IsSmallValue, so it is non-tiny here.)
		using Real = typename Eigen::NumTraits<T>::Real;
		const Real d = Eigen::NumTraits<T>::dummy_precision();
		return Eigen::numext::abs2(numerator) * (d * d) >= Eigen::numext::abs2(denomenator);
	}

	/// \brief Outcome of a matrix-decomposition sanity check (e.g. LU pivot magnitude).
	enum class MatrixSuccessCode
	{
		Success,      ///< The decomposition looks healthy.
		LargeChange,  ///< A pivot changed by a suspiciously large amount.
		SmallValue    ///< A pivot is suspiciously small (near-singular).
	};

	/**
	\brief Check the diagonal elements of an LU decomposition for small values and large ratios.

	\return Success if things are ok.  LargeChange or SmallValue if one is found.

	This function requires a square non-empty matrix.

	\tparam Derived Matrix type from Eigen.
	*/
	template <typename Derived>
	MatrixSuccessCode LUPartialPivotDecompositionSuccessful(Eigen::MatrixBase<Derived> const& LU)
	{
		#ifndef BERTINI_DISABLE_ASSERTS
			assert(LU.rows()==LU.cols() && "non-square matrix in LUPartialPivotDecompositionSuccessful");
			assert(LU.rows()>0 && "empty matrix in LUPartialPivotDecompositionSuccessful");
		#endif

			// this loop won't test entry (0,0).  it's tested separately after.
		for (Eigen::Index ii = LU.rows()-1; ii > 0; ii--)
		{
			if (IsSmallValue(LU(ii,ii)))
			{
				return MatrixSuccessCode::SmallValue;
			}

			if (IsLargeChange(LU(ii-1,ii-1),LU(ii,ii)))
			{
				return MatrixSuccessCode::LargeChange;
			}
		}

		// this line is the reason for the above assert on non-empty matrix.
		if (IsSmallValue(LU(0,0)))
		{
			return MatrixSuccessCode::SmallValue;
		}

		return MatrixSuccessCode::Success;
	}

	/**
	\brief Make a Kahan matrix with a given number type.
	*/
	template <typename NumberType>
	Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic> KahanMatrix(unsigned int mat_size, NumberType c)
	{
		using std::sqrt;
		NumberType s, scale(1.0);
		s = sqrt( (NumberType(1.0)-c) * (NumberType(1.0)+c) );


		Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic> A(mat_size,mat_size);


		for (unsigned int ii=0; ii<mat_size; ii++) {
			for (unsigned int jj=0; jj<ii; jj++) {
				A(ii,jj) = NumberType(0.0);
			}
			for (unsigned int jj=ii; jj<mat_size; jj++) {
				A(ii,jj) = -c * NumberType(1.0);
			}
		}


		for (unsigned int ii=0; ii<mat_size; ii++) {
			A(ii,ii) += NumberType(1)+c;
		}

		for (unsigned int jj=0; jj<mat_size; jj++){
			for (unsigned int kk=0;kk<mat_size;kk++){
				A(jj,kk) *= scale;
			}
			scale *= s;
		}

		for (unsigned int jj=0;jj<mat_size;jj++){
			for (unsigned int kk=0;kk<mat_size;kk++){
				A(kk,jj)/= NumberType(jj) + NumberType(1);
			}
		}
		return A;
	}


	/**
	\brief Make a random matrix of units (numbers with norm 1).

	\return The random matrix of units.

	\param rows The number of rows.
	\param cols The number of columns.
	\tparam NumberType the type of number to fill the matrix with.
	*/
	template <typename NumberType>
	inline
	Mat<NumberType> RandomOfUnits(unsigned int rows, unsigned int cols)
	{
		return Mat<NumberType>(rows,cols).unaryExpr([](NumberType const& /*x*/) { return RandomUnit<NumberType>(); });
	}

	/**
	\brief Make a random vector of units (numbers with norm 1).

	\return The random vector of units.

	\param size The length of the vector.
	\tparam NumberType the type of number to fill the vector with.
	*/
	template <typename NumberType>
	inline
	Vec<NumberType> RandomOfUnits(unsigned int size)
	{
		return Vec<NumberType>(size).unaryExpr([](NumberType const&) { return RandomUnit<NumberType>(); });
	}

}


// // the following code comes from
// // https://stackoverflow.com/questions/18382457/eigen-and-boostserialize
// // and adds support for serialization of Eigen types
// //
// // question asked by user Gabriel and answered by Shmuel Levine
// // answer code repeated here verbatim.
// // please update this comment if this code is changed,
// // and post the modifications to the above referenced post on SO.
// namespace boost{
//     namespace serialization{

//         template<   class Archive,
//                     class S,
//                     int Rows_,
//                     int Cols_,
//                     int Ops_,
//                     int MaxRows_,
//                     int MaxCols_>
//         inline void save(
//             Archive & ar,
//             const Eigen::Matrix<S, Rows_, Cols_, Ops_, MaxRows_, MaxCols_> & g,
//             const unsigned int version)
//             {
//                 int rows = g.rows();
//                 int cols = g.cols();

//                 ar & rows;
//                 ar & cols;
//                 ar & boost::serialization::make_array(g.data(), rows * cols);
//             }

//         template<   class Archive,
//                     class S,
//                     int Rows_,
//                     int Cols_,
//                     int Ops_,
//                     int MaxRows_,
//                     int MaxCols_>
//         inline void load(
//             Archive & ar,
//             Eigen::Matrix<S, Rows_, Cols_, Ops_, MaxRows_, MaxCols_> & g,
//             const unsigned int version)
//         {
//             int rows, cols;
//             ar & rows;
//             ar & cols;
//             g.resize(rows, cols);
//             ar & boost::serialization::make_array(g.data(), rows * cols);
//         }

//         template<   class Archive,
//                     class S,
//                     int Rows_,
//                     int Cols_,
//                     int Ops_,
//                     int MaxRows_,
//                     int MaxCols_>
//         inline void serialize(
//             Archive & ar,
//             Eigen::Matrix<S, Rows_, Cols_, Ops_, MaxRows_, MaxCols_> & g,
//             const unsigned int version)
//         {
//             split_free(ar, g, version);
//         }


//     } // namespace serialization
// } // namespace boost




#endif
