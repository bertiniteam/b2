//This file is part of Bertini 2.
//
// bertini2/random.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
// bertini2/random.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  bertini2/random.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file  bertini2/random.hpp 

\brief stuff for generating random numbers
*/

#ifndef BERTINI_RANDOM_HPP
#define BERTINI_RANDOM_HPP




#include "bertini2/mpfr_complex.hpp"
#include "bertini2/records/derive.hpp"
#include <boost/random.hpp>
#include <random>
#include <cstdint>


namespace bertini
{

	/**
	Returns this thread's canonical mt19937 engine.  Every random draw of every
	type routes through here — integer/rational (RandomInt/RandomRat), double
	(rand_complex/RandReal), and multiprecision (RandomMp and everything built on
	it) — so SetGlobalSeed()/ReseedThisThread() control them all uniformly.
	*/
	std::mt19937& ThreadEngine();

	/**
	Set the global RNG seed.  seed == 0 draws from std::random_device and stores
	the effective (non-zero) seed so it can be retrieved and reproduced later.
	Must be called before system construction (gamma, patch, TD-constants) to make
	those setup draws deterministic.
	*/
	void SetGlobalSeed(unsigned long seed);

	/**
	Returns the effective global seed.  If SetGlobalSeed has never been called, draws
	from entropy on first call and caches the result.
	*/
	unsigned long GetGlobalSeed();

	/**
	Reseed this thread's engine deterministically from the global seed mixed with
	stream_key.  Call at the top of each TrackSinglePath* with soln_ind as the key
	so per-step random draws (condition-number probe, PSEG rand-vector) are
	path-indexed and mode-independent.
	*/
	void ReseedThisThread(uint64_t stream_key);

	/**
	Derive a distinct, deterministic child seed for a worker rank from the global (master) seed.
	The MPI manager computes one per worker and hands it over; the worker calls SetGlobalSeed(child),
	giving every process its own non-overlapping deterministic stream -- all reproducible from the one
	user seed, and no two processes ever generate the same random value.
	*/
	unsigned long DerivedWorkerSeed(uint64_t worker_index);


	/**
	Generate a random integer number between -10^digits and 10^digits
	*/
	template <unsigned long digits = 50>
	inline
	mpz_int RandomInt()
	{
		// pinned draw (b2rand/1, ADR-0044): uniform on [-2^(digits/log10(2)), +same]
		return records::ThreadDrawStream().IntSymmetric(mpz_int(1) << digits*1000L/301L);
	}


	/**
	Generate a random rational number with numerator and denomenator between -10^digits and 10^digits
	*/
	template <unsigned long digits = 50>
	mpq_rational RandomRat()
	{
		// pinned draws (b2rand/1, ADR-0044).  A zero denominator is redrawn (deterministically):
		// the legacy draw could in principle hand mpq a denominator of 0.
		mpz_int const bound = mpz_int(1) << digits*1000L/301L;
		auto& stream = records::ThreadDrawStream();
		mpz_int const num = stream.IntSymmetric(bound);
		mpz_int den = stream.IntSymmetric(bound);
		while (den == 0)
			den = stream.IntSymmetric(bound);
		return mpq_rational(num, den);
	}


	/**
	 Produce a random number with at length_in_digits non-zero digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 */
	template <unsigned int length_in_digits>
	real_mp RandomMp()
	{
		// Pinned draw (b2rand/1, ADR-0044) from the single per-thread stream shared by
		// every random type, so SetGlobalSeed()/ReseedThisThread() control them all
		// uniformly.  Uniform on [0,1) at length_in_digits digits, bit-identical on
		// every platform.
		real_mp a{records::ThreadDrawStream().UnitRealMp(length_in_digits)};
		return a;
	}
	
	/**
	 a templated function for producing random numbers in the unit interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 \param a the number which will be assigned in this call
	 */
	template <unsigned int length_in_digits>
	void RandomMpAssign(real_mp & a)
	{	
		a = RandomMp<length_in_digits>();
	}

	/**
	 a templated function for producing random numbers in a specified interval, of a given number of digits.
	 
	 \tparam length_in_digits The length of the desired random number
	 
	 \param left The left bound.
	 \param right The right bound.
	 */
	template <unsigned int length_in_digits>
	real_mp RandomMp(const real_mp & left, const real_mp & right)
	{
		return (right-left)*RandomMp<length_in_digits>()+left;
	}



	/**
	 \brief create a random number, at the current default precision
	 */
	real_mp RandomMp();

	/**
	 \brief create a random number, at the specified precision

	 \param num_digits the precision that you desire.  
	 */
	real_mp RandomMp(unsigned num_digits);

	/**
	 \brief create a random number in a given interval, at the current default precision
	*/
	real_mp RandomMp(const real_mp & a, const real_mp & b);

	/**
	 \brief create a random number in a given interval, at the specified precision
	*/
	real_mp RandomMp(const real_mp & a, const real_mp & b, unsigned num_digits);

	/**
	 \brief Set an existing real_mp to a random number, to a given precision.  

	 This function is how to get random numbers at a precision different from the current default.
	 */
	void RandomMpAssign(real_mp & a, unsigned num_digits);

	

	
} // re: namespace bertini




namespace bertini{

namespace multiprecision{


using complex = bertini::complex_mp;  ///< Shorthand for the multiprecision complex type within this namespace.
using bertini::RandomMp;



	/**
	 Assign to a random real number \f$\in [-1,\,1]\f$, to current default precision. 
	 */
	inline 
	void RandomRealAssign(complex & a, unsigned num_digits)
	{
		auto cached = ThreadPrecision();
		SetThreadPrecision(num_digits);
		complex temp(RandomMp(real_mp(-1),real_mp(1),num_digits)); // ,0
		a.swap(temp);
		SetThreadPrecision(cached);
	}

	/**
	 Produce a random real number \f$\in [-1,\,1]\f$, to current default precision. 
	 */
	inline complex RandomReal()
	{
		return complex(RandomMp(real_mp(-1),real_mp(1))); // ,0
	}
	
	/**
	 Produce a random real number \f$\in [-1,\,1]\f$, to specified precision. 
	 */
	inline complex RandomReal(unsigned num_digits)
	{
		auto cached = ThreadPrecision();
		SetThreadPrecision(num_digits);
		auto result = complex(RandomMp(real_mp(-1),real_mp(1),num_digits));// ,0
		SetThreadPrecision(cached);
		return result;
	}



	
		


	/**
	 Produce a random complex number, to default precision.
	 */
	inline complex rand()
	{
		return complex( RandomMp(real_mp(-1),real_mp(1)), RandomMp(real_mp(-1),real_mp(1)) );
	}


	/**
	 Produce a random unit complex number, to default precision.
	 */
	inline complex rand_unit()
	{
		complex returnme( RandomMp(real_mp(-1),real_mp(1)), RandomMp(real_mp(-1),real_mp(1)) );
		return returnme / abs(returnme);   // normalize to modulus 1 (NOT sqrt(abs), which left modulus sqrt|z|)
	}

	/// \brief Produce a random unit-modulus complex number, to default precision.
	inline complex RandomUnit()
	{
		return rand_unit();
	}


	/**
	 Produce a random complex number whose modulus is pulled toward 1, to default precision.

	 Draw a box-uniform complex z (real, imag each in [-1,1]) and divide by sqrt(|z|), so the
	 result has modulus sqrt(|z|): bounded away from both 0 and infinity, but NOT collapsed onto
	 the unit circle (that would be z/abs(z) -- see rand_unit).  This is how Bertini 1 generates
	 linear-form coefficients, and it avoids the heavy-tailed scaling of a ratio-of-integers draw
	 (RandomRat: numerator/denominator each uniform, so the modulus has fat log-tails).
	 */
	inline complex rand_bounded_modulus()
	{
		complex z( RandomMp(real_mp(-1),real_mp(1)), RandomMp(real_mp(-1),real_mp(1)) );
		auto m = abs(z);
		while (m == 0)   // measure-zero, but a zero coefficient is degenerate -- redraw
		{
			z = complex( RandomMp(real_mp(-1),real_mp(1)), RandomMp(real_mp(-1),real_mp(1)) );
			m = abs(z);
		}
		return z / sqrt(m);
	}

	/// \brief Produce a random complex number whose modulus is pulled toward 1 (away from 0 and infinity), at default precision.
	inline complex RandomComplexBoundedModulus()
	{
		return rand_bounded_modulus();
	}

	/// \brief Assign to a a random complex number whose modulus is pulled toward 1, at the given precision.
	inline
	void RandomComplexBoundedModulusAssign(complex & a, unsigned num_digits)
	{
		auto cached = ThreadPrecision();
		SetThreadPrecision(num_digits);
		a.precision(num_digits);

		complex z( RandomMp(real_mp(-1),real_mp(1),num_digits), RandomMp(real_mp(-1),real_mp(1),num_digits) );
		auto m = abs(z);
		while (m == 0)
		{
			z = complex( RandomMp(real_mp(-1),real_mp(1),num_digits), RandomMp(real_mp(-1),real_mp(1),num_digits) );
			m = abs(z);
		}
		a = std::move(z / sqrt(m));
		SetThreadPrecision(cached);
	}

	/// \brief Produce a random complex number whose modulus is pulled toward 1, at the given precision.
	inline
	complex RandomComplexBoundedModulus(unsigned num_digits)
	{
		complex a;
		RandomComplexBoundedModulusAssign(a, num_digits);
		return a;
	}


	/**
	 Produce a random REAL number whose modulus is pulled toward 1, to default precision.

	 The real-line analog of rand_bounded_modulus: draw x box-uniform in [-1,1] (imaginary part 0)
	 and divide by sqrt(|x|), giving sign(x)*sqrt(|x|) -- bounded away from both 0 and infinity, the
	 same Bertini-1 recipe the (complex) start-system and patch coefficients use, but kept REAL so a
	 real patch does not complexify a real path.
	 */
	inline complex rand_real_bounded_modulus()
	{
		complex z( RandomMp(real_mp(-1),real_mp(1)) );   // imaginary part 0
		auto m = abs(z);
		while (m == 0)   // measure-zero, but a zero coefficient is degenerate -- redraw
		{
			z = complex( RandomMp(real_mp(-1),real_mp(1)) );
			m = abs(z);
		}
		return z / sqrt(m);   // stays real (imag 0 / real = 0)
	}

	/// \brief Produce a random REAL number (imaginary part 0) whose modulus is pulled toward 1, at default precision.
	inline complex RandomRealBoundedModulus()
	{
		return rand_real_bounded_modulus();
	}

	/// \brief Assign to a a random REAL number (imaginary part 0) whose modulus is pulled toward 1, at the given precision.
	inline
	void RandomRealBoundedModulusAssign(complex & a, unsigned num_digits)
	{
		auto cached = ThreadPrecision();
		SetThreadPrecision(num_digits);
		a.precision(num_digits);

		complex z( RandomMp(real_mp(-1),real_mp(1),num_digits) );   // imaginary part 0
		auto m = abs(z);
		while (m == 0)
		{
			z = complex( RandomMp(real_mp(-1),real_mp(1),num_digits) );
			m = abs(z);
		}
		a = std::move(z / sqrt(m));   // stays real
		SetThreadPrecision(cached);
	}

	/// \brief Produce a random REAL number (imaginary part 0) whose modulus is pulled toward 1, at the given precision.
	inline
	complex RandomRealBoundedModulus(unsigned num_digits)
	{
		complex a;
		RandomRealBoundedModulusAssign(a, num_digits);
		return a;
	}

	/// \brief Assign a random complex number to a, at the given precision.
	inline
	void rand_assign(complex & a, unsigned num_digits)
	{
		auto cached = ThreadPrecision();
		SetThreadPrecision(num_digits);
		
		complex_mp temp( RandomMp(num_digits), RandomMp(num_digits) );
		a = std::move(temp);
		SetThreadPrecision(cached);
	}

	/// \brief Assign a random complex number to a, at the given precision.
	inline
	void RandomComplexAssign(complex & a, unsigned num_digits)
	{
		rand_assign(a,num_digits);
	}

	/// \brief Produce a random complex number, at the given precision.
	inline
	complex RandomComplex(unsigned num_digits)
	{
		complex z;
		RandomComplexAssign(z, num_digits);
		return z;
	}

	/// \brief Assign a random unit-modulus complex number to a, at the given precision.
	inline
	void RandomUnitAssign(complex & a, unsigned num_digits)
	{
		auto cached = ThreadPrecision();
		SetThreadPrecision(num_digits);
		a.precision(num_digits);
		
		complex temp(RandomMp(num_digits),RandomMp(num_digits));
		a = std::move(temp/abs(temp));   // normalize to modulus 1 (NOT sqrt(abs))
		SetThreadPrecision(cached);
	}

	/// \brief Produce a random unit-modulus complex number, at the given precision.
	inline
	complex RandomUnit(unsigned num_digits)
	{
		complex a;
		RandomUnitAssign(a,num_digits);
		return a;
	}

}  // namespace multiprecision


} // namespaces
// }



#endif




