//This file is part of Bertini 2.
//
//random.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//random.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with random.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file random.cpp 

\brief stuff to make random numbers in bertini2
*/

#include "bertini2/random.hpp"

#include <random>
#include <atomic>
#include <cstdint>

namespace bertini {

namespace {

inline uint64_t splitmix64(uint64_t x)
{
	x += 0x9e3779b97f4a7c15ULL;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
	return x ^ (x >> 31);
}

std::atomic<unsigned long> g_global_seed{0};

// Domain tags for stream derivation.  Distinct domains (and distinct indices within a domain)
// produce distinct engine states, so no two streams ever coincide and no path/worker stream can
// reproduce the setup stream.  (Setup = the stream that generates gamma / start coefficients / patch.)
constexpr uint64_t kDomainSetup  = 0x5e7400000000ULL; // "setup"  -- SetGlobalSeed
constexpr uint64_t kDomainPath   = 0x9a7400000000ULL; // "path"   -- per-path tracking streams
constexpr uint64_t kDomainWorker = 0x107ce000000ULL;  // "worker" -- per-rank child seeds (MPI)

// Seed the engine from the FULL (master, domain, index) tuple via std::seed_seq.  Using the whole
// 64-bit words (not a uint32_t truncation) means distinct tuples set distinct mt19937 states, so
// streams cannot collide modulo 2^32 (the old bug).
inline void SeedEngine(std::mt19937& eng, uint64_t master, uint64_t domain, uint64_t index)
{
	std::seed_seq seq{
		static_cast<uint32_t>(master),        static_cast<uint32_t>(master >> 32),
		static_cast<uint32_t>(domain),        static_cast<uint32_t>(domain >> 32),
		static_cast<uint32_t>(index),         static_cast<uint32_t>(index >> 32)
	};
	eng.seed(seq);
}

} // anon namespace

// one mt19937 per thread — seeded from entropy on first use (overwritten deterministically by
// SetGlobalSeed / ReseedThisThread before any draw in a seeded run).
thread_local std::mt19937 g_thread_engine{std::random_device{}()};

std::mt19937& ThreadEngine() { return g_thread_engine; }

unsigned long GetGlobalSeed()
{
	unsigned long s = g_global_seed.load(std::memory_order_relaxed);
	if (s == 0) {
		std::random_device rd;
		s = static_cast<unsigned long>(rd());
		if (s == 0) s = 1;
		unsigned long expected = 0;
		if (!g_global_seed.compare_exchange_strong(expected, s, std::memory_order_relaxed))
			s = g_global_seed.load(std::memory_order_relaxed);
	}
	return s;
}

void SetGlobalSeed(unsigned long seed)
{
	if (seed == 0) {
		std::random_device rd;
		seed = static_cast<unsigned long>(rd());
		if (seed == 0) seed = 1;
	}
	g_global_seed.store(seed, std::memory_order_relaxed);
	// the setup stream: domain = setup, index = 0.  Path/worker streams use other domains, so none
	// of them can ever reproduce this stream (the old ReseedThisThread(0) == SetGlobalSeed collision).
	SeedEngine(g_thread_engine, static_cast<uint64_t>(seed), kDomainSetup, 0);
}

void ReseedThisThread(uint64_t stream_key)
{
	// per-path / per-thread stream: domain = path, index = stream_key.  Deterministic from the global
	// seed and distinct for every stream_key (and distinct from the setup stream).
	SeedEngine(g_thread_engine, static_cast<uint64_t>(GetGlobalSeed()), kDomainPath, stream_key);
}

// Derive a distinct, deterministic child seed for a worker rank, from the master seed.  The manager
// computes these and hands one to each worker, which then calls SetGlobalSeed(child) -- so every
// process has its own non-overlapping deterministic stream, all reproducible from the one user seed.
unsigned long DerivedWorkerSeed(uint64_t worker_index)
{
	uint64_t s = static_cast<uint64_t>(GetGlobalSeed());
	uint64_t h = splitmix64(s ^ kDomainWorker ^ splitmix64(worker_index));
	if (h == 0) h = 1; // SetGlobalSeed treats 0 as "draw from entropy"; avoid that
	return static_cast<unsigned long>(h);
}



	mpfr_float RandomMp()
	{
		// ThreadPrecision (thread-local) rather than DefaultPrecision (global):
		// random numbers are generated during tracking, which may run on a
		// std::thread worker whose precision differs from the global default.
		return RandomMp(bertini::ThreadPrecision());
	}

	mpfr_float RandomMp(unsigned num_digits)
	{
		
		mpfr_float a;
		if (num_digits<=50)
			a = std::move(RandomMp<50>());
		else if (num_digits<=100)
			a = std::move(RandomMp<100>());
		else if (num_digits<=200)
			a = std::move(RandomMp<200>());
		else if (num_digits<=400)
			a = std::move(RandomMp<400>());
		else if (num_digits<=800)
			a = std::move(RandomMp<800>());
		else if (num_digits<=1600)
			a = std::move(RandomMp<1600>());
		else if (num_digits<=3200)
			a = std::move(RandomMp<3200>());
		else if (num_digits<=6400)
			a = std::move(RandomMp<6400>());
		else if (num_digits<=8000)
			a = std::move(RandomMp<8000>());
		else if (num_digits<=10000)
			a = std::move(RandomMp<10000>());
		else if (num_digits<=12000)
			a = std::move(RandomMp<12000>());
		else if (num_digits<=14000)
			a = std::move(RandomMp<14000>());
		else if (num_digits<=16000)
			a = std::move(RandomMp<16000>());
		else if (num_digits<=18000)
			a = std::move(RandomMp<18000>());
		else if (num_digits<=20000)
			a = std::move(RandomMp<20000>());
		else if (num_digits<=40000)
			a = std::move(RandomMp<40000>());
		else
			throw std::out_of_range("requesting random long number of digits -- higher than 40000.  this throw can be remedied by adding more cases to the generating function RandomMp in random.cpp.  If you have a better solution to this problem, please write the authors of this software.");
		a.precision(num_digits);
		return a;
	}


	void RandomMpAssign(mpfr_float & a, unsigned num_digits)
	{
		

		mpfr_float temp;
		temp = RandomMp(num_digits);
		a = std::move(temp);
	}




	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b)
	{
		// see RandomMp() above for why ThreadPrecision rather than DefaultPrecision
		return RandomMp(a,b,bertini::ThreadPrecision());
	}

	mpfr_float RandomMp(const mpfr_float & a, const mpfr_float & b, unsigned num_digits)
	{
		
		mpfr_float result;
		if (num_digits<=50)
			result = std::move(RandomMp<50>(a,b));
		else if (num_digits<=100)
			result = std::move(RandomMp<100>(a,b));
		else if (num_digits<=200)
			result = std::move(RandomMp<200>(a,b));
		else if (num_digits<=400)
			result = std::move(RandomMp<400>(a,b));
		else if (num_digits<=800)
			result = std::move(RandomMp<800>(a,b));
		else if (num_digits<=1600)
			result = std::move(RandomMp<1600>(a,b));
		else if (num_digits<=3200)
			result = std::move(RandomMp<3200>(a,b));
		else if (num_digits<=6400)
			result = std::move(RandomMp<6400>(a,b));
		else if (num_digits<=8000)
			result = std::move(RandomMp<8000>(a,b));
		else if (num_digits<=10000)
			result = std::move(RandomMp<10000>(a,b));
		else if (num_digits<=12000)
			result = std::move(RandomMp<12000>(a,b));
		else if (num_digits<=14000)
			result = std::move(RandomMp<14000>(a,b));
		else if (num_digits<=16000)
			result = std::move(RandomMp<16000>(a,b));
		else if (num_digits<=18000)
			result = std::move(RandomMp<18000>(a,b));
		else if (num_digits<=20000)
			result = std::move(RandomMp<20000>(a,b));
		else if (num_digits<=40000)
			result = std::move(RandomMp<40000>(a,b));
		else
			throw std::out_of_range("requesting random long number of digits -- higher than 40000.  this throw can be remedied by adding more cases to the generating function RandomMp in random.cpp.  If you have a better solution to this problem, please write the authors of this software.");
		result.precision(num_digits);
		return result;
	}






	


} // namespace bertini
