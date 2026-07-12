//This file is part of Bertini 2.
//
//num_traits.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//num_traits.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with num_traits.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file num_traits.hpp 

\brief Provides an Eigen-like NumTraits struct for querying traits of a number type.  

The bertini::NumTraits struct provides NumDigits and NumFuzzyDigits functions.
*/

#ifndef BERTINI_NUM_TRAITS_HPP
#define BERTINI_NUM_TRAITS_HPP

#include <random>
#include "bertini2/records/draw_functions.hpp"
#include <complex>
#include <cmath>
#include "bertini2/mpfr_complex.hpp"
#include "bertini2/random.hpp"


namespace bertini
{
	/// \brief The numeric type used to hold error/tolerance magnitudes throughout Bertini2 --
	/// tracking tolerances, Newton residuals, condition numbers, same-point tolerances.  Declared
	/// here in the foundational num_traits header (included nearly everywhere) so tolerance-typed
	/// parameters can spell this alias instead of a bare double without any header adding a new
	/// include -- keeping recompilation frequency and scope unchanged.
	using NumErrorT = double;

	/// \brief Get a random complex number of unit modulus, in number type T.
	template<typename T>
	T RandomUnit();

	/// \brief Numeric traits for a bertini number type: digit counts, conversions, and companion types.
	template<typename T>
	struct NumTraits
	{};



	/// \brief Numeric traits for the built-in double type.
	template <> struct NumTraits<double>
	{
		/// \brief The number of significant digits (16 for double).
		inline static unsigned NumDigits()
		{
			return 16;
		}

		/// \brief The number of digits to trust for fuzzy comparisons.
		inline static unsigned NumFuzzyDigits()
		{
			return 14;
		}

		/// \brief Convert a tracking tolerance to a number of significant digits.
		inline
		static unsigned TolToDigits(NumErrorT tol)
		{
			return static_cast<unsigned>(ceil(-log10(tol)));
		}

		/// \brief Parse a double from a string.
		inline static
		double FromString(std::string const& s)
		{
			return boost::lexical_cast<double>(s);
		}

		/// \brief Convert an exact rational to a double (precision ignored).
		inline static
		double FromRational(mpq_rational const& n, unsigned /* precision */)
		{
			return double(n);
		}

		using Real = double;  ///< The real companion type.
		using Complex = complex_dbl;  ///< The complex companion type.
	};


	/// \brief Numeric traits for double-precision complex numbers.
	template <> struct NumTraits<complex_dbl >
	{
		/// \brief The number of significant digits (16).
		inline static unsigned NumDigits()
		{
			return 16;
		}

		/// \brief The number of digits to trust for fuzzy comparisons.
		inline static unsigned NumFuzzyDigits()
		{
			return 14;
		}

		/// \brief Parse a complex number from a single string.
		inline static
		complex_dbl FromString(std::string const& s)
		{
			return boost::lexical_cast<complex_dbl>(s);
		}

		/// \brief Build a complex number from real and imaginary string parts.
		inline static
		complex_dbl FromString(std::string const& s, std::string const& t)
		{
			return complex_dbl(boost::lexical_cast<double>(s),boost::lexical_cast<double>(t));
		}

		/// \brief Convert an exact rational to a complex double (precision ignored).
		inline static
		complex_dbl FromRational(mpq_rational const& n, unsigned /* precision */)
		{
			return complex_dbl(static_cast<double>(n),0);
		}

		using Real = double;  ///< The real companion type.
		using Complex = complex_dbl;  ///< The complex companion type.
	};


	/// \brief The number of digits by which adaptive precision is increased at each step.
	inline
	unsigned PrecisionIncrement()
	{
		return 10;
	}

	/// \brief The number of digits considered "double precision" (16).
	constexpr
	unsigned DoublePrecision()
	{
		return 16;
	}

	/// \brief The lowest multiple-precision digit count used (20).
	constexpr
	unsigned LowestMultiplePrecision()
	{
		return 20;
	}

	/// \brief The maximum precision (digits) the library will escalate to.
	constexpr
	unsigned MaxPrecisionAllowed()
	{
		return 1000;
	}
	
	/**
	\brief Get the precision of a number.

	For doubles, this is trivially 16.
	*/
	inline
	unsigned Precision(double)
	{
		return DoublePrecision();
	}

	/**
	For complex doubles, throw if the requested precision is not DoublePrecision.
	*/
	inline
	void Precision(double, unsigned prec)
	{
		if (prec!=DoublePrecision())
		{
			std::stringstream err_msg;
			err_msg << "trying to change precision of a double to " << prec;
			throw std::runtime_error(err_msg.str());
		}
	}

	/**
	\brief Get the precision of a number.

	For complex doubles, this is trivially 16.
	*/
	inline
	unsigned Precision(complex_dbl)
	{
		return DoublePrecision();
	}

	/**
	For complex doubles, throw if the requested precision is not DoublePrecision.
	*/
	inline
	void Precision(complex_dbl, unsigned prec)
	{
		if (prec!=DoublePrecision())
		{
			std::stringstream err_msg;
			err_msg << "trying to change precision of a double to " << prec;
			throw std::runtime_error(err_msg.str());
		}
	}

	/// \brief Get a random double-precision complex number (not unit modulus; see RandomUnit).
	inline
	complex_dbl rand_complex()
	{
		using std::abs;
		using std::sqrt;
		// DRAW ORDER IS A CONTRACT (cross-platform reproducibility, 2026-07-03): the two
		// component draws are sequenced EXPLICITLY -- real first, then imaginary.  Never
		// put two draws in one full-expression: C++ argument evaluation order is
		// unspecified, and gcc really does order them differently on x86_64 vs aarch64
		// (found as an architecture-split seeded-homotopy digest: every (re, im) pair of
		// every seeded coefficient was transposed between the two).
		double const re = records::DrawSymmetricDouble();
		double const im = records::DrawSymmetricDouble();
		complex_dbl returnme(re, im);
		return returnme / sqrt( abs(returnme));
	}

	/// \brief Get a random double-precision complex number of unit modulus.
	template <> inline
	complex_dbl RandomUnit<complex_dbl >()
	{
		// draw order is a contract: real first, then imaginary (see rand_complex)
		double const re = records::DrawSymmetricDouble();
		double const im = records::DrawSymmetricDouble();
		complex_dbl returnme(re, im);
		return returnme / abs(returnme);
	}

	/// \brief Get a random multiprecision complex number of unit modulus.
	template <>
	inline
	complex_mp RandomUnit<complex_mp>()
	{
		return multiprecision::RandomUnit();
	}

	/// \brief Get a random multiprecision real number of unit modulus (i.e. +1 or -1).
	/// The real analog of the complex unit draw, so RandomConjugateOrthonormalMatrix<real_mp>
	/// (QR of a matrix of real units) yields a real orthogonal matrix.
	template <>
	inline
	real_mp RandomUnit<real_mp>()
	{
		return real_mp( records::DrawSymmetricDouble() < 0.0 ? -1 : 1 );
	}
}// re: namespace bertini












namespace bertini {

	
	
	/// \brief Numeric traits for multiprecision reals (digit counts follow the thread precision).
	template <> struct NumTraits<real_mp>
	{
		/// \brief The number of significant digits (the current thread precision).
		inline static unsigned NumDigits()
		{
			return ThreadPrecision();
		}

		/// \brief The number of digits to trust for fuzzy comparisons.
		inline static unsigned NumFuzzyDigits()
		{
			return ThreadPrecision()-3;
		}

		/// \brief Convert a tracking tolerance to a number of significant digits.
		inline
		static unsigned TolToDigits(real_mp tol)
		{
			real_mp b = ceil(-log10(tol));
			return b.convert_to<unsigned int>();
		}

		/// \brief Parse a multiprecision real from a string.
		inline static
		real_mp FromString(std::string const& s)
		{
			return real_mp(s);
		}

		/// \brief Convert an exact rational to a multiprecision real at the given precision.
		inline static
		real_mp FromRational(mpq_rational const& n, unsigned precision)
		{
			return real_mp(n,precision);
		}

		using Real = real_mp;  ///< The real companion type.
		using Complex = complex_mp;  ///< The complex companion type.
	};



	/// \brief Numeric traits for multiprecision complex numbers (digit counts follow the thread precision).
	template <> struct NumTraits<complex_mp>
	{
		/// \brief The number of significant digits (the current thread precision).
		inline static unsigned NumDigits()
		{
			return ThreadPrecision();
		}

		/// \brief Parse a complex number from a single string.
		inline static
		complex_mp FromString(std::string const& s)
		{
			return complex_mp(s);
		}

		/// \brief Build a complex number from real and imaginary string parts.
		inline static
		complex_mp FromString(std::string const& s, std::string const& t)
		{
			return complex_mp(s,t);
		}

		/// \brief Convert an exact rational to a multiprecision complex number at the given precision.
		inline static
		complex_mp FromRational(mpq_rational const& n, unsigned precision)
		{
			return complex_mp(n,0,precision);
		}

		using Real = real_mp;  ///< The real companion type.
		using Complex = complex_mp;  ///< The complex companion type.
	};

	/// \brief Numeric traits for exact rationals (provides decimal-string parsing).
	template <> struct NumTraits<mpq_rational>
	{
		/**
		\brief Parse a decimal string to an exact mpq_rational.

		Handles optional sign, decimal point, and scientific notation (e/E).
		Unlike the mpq_rational(string) constructor (which expects "p/q" or integer
		format), this accepts decimal strings like "0.647" → 647/1000 exactly.

		GMP treats a leading '0' as an octal prefix when base=0, so leading zeros
		are stripped before mpz_int construction.  "0.8" → digits "08" → "8".
		*/
		inline static
		mpq_rational FromString(std::string const& str)
		{
			std::string s = str;
			bool negative = false;
			if (!s.empty() && s[0] == '-') { negative = true; s = s.substr(1); }
			else if (!s.empty() && s[0] == '+') { s = s.substr(1); }

			int exp_shift = 0;
			auto e_pos = s.find_first_of("eE");
			if (e_pos != std::string::npos) {
				exp_shift = std::stoi(s.substr(e_pos + 1));
				s = s.substr(0, e_pos);
			}

			auto dot_pos = s.find('.');
			int decimal_places = 0;
			if (dot_pos != std::string::npos) {
				decimal_places = static_cast<int>(s.size()) - static_cast<int>(dot_pos) - 1;
				s.erase(dot_pos, 1);
			}

			if (s.empty() || s.find_first_not_of('0') == std::string::npos) {
				s = "0";
			} else {
				s = s.substr(s.find_first_not_of('0'));
			}
			mpz_int numer(s);
			if (negative) numer = -numer;

			int net_exp = decimal_places - exp_shift;
			if (net_exp > 0) {
				mpz_int denom = 1;
				for (int i = 0; i < net_exp; ++i) denom *= 10;
				return mpq_rational(numer, denom);
			} else if (net_exp < 0) {
				mpz_int mult = 1;
				for (int i = 0; i < -net_exp; ++i) mult *= 10;
				return mpq_rational(numer * mult, 1);
			} else {
				return mpq_rational(numer, 1);
			}
		}
	};

}

#endif


