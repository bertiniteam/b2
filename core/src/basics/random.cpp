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

/**
\file random.cpp

\brief stuff to make random numbers in bertini2
*/

#include "bertini2/random.hpp"
#include "bertini2/records/derive.hpp"

#include <random>
#include <atomic>
#include <cstdint>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace bertini {

namespace {

inline uint64_t splitmix64(uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// The global seed, and whether one has been chosen.  Every value is a seed, 0 included, so
// "not chosen yet" is its own flag; the mutex only guards the one-time entropy draw.
std::atomic<unsigned long> g_global_seed{0};
std::atomic<bool> g_global_seed_chosen{false};
std::mutex g_global_seed_mutex;

unsigned long EntropySeed()
{
    std::random_device rd;
    return static_cast<unsigned long>(rd());
}

// Domain tags for stream derivation.  Distinct domains (and distinct indices within a domain)
// produce distinct engine states, so no two streams ever coincide and no path/worker stream can
// reproduce the setup stream.  (Setup = the stream that generates gamma / start coefficients / patch.)
constexpr uint64_t kDomainSetup  = 0x5e7400000000ULL; // "setup"  -- SetGlobalSeed
constexpr uint64_t kDomainPath   = 0x9a7400000000ULL; // "path"   -- per-path tracking streams
constexpr uint64_t kDomainWorker = 0x107ce000000ULL;  // "worker" -- per-rank child seeds (MPI)
constexpr uint64_t kDomainSolve  = 0x501e000000000ULL; // "solve"  -- effective per-solve seeds

// Ordinal the next AdvanceSolveSeed call takes under the current master; reset by SetGlobalSeed
// so (master, ordinal) is well-defined.  Since a solve rekeys to its own effective seed,
// consecutive seedless solves form a deterministic seed CHAIN from the initial master.
std::atomic<uint64_t> g_solve_ordinal{0};

} // anon namespace

unsigned long GetGlobalSeed()
{
    if (g_global_seed_chosen.load(std::memory_order_acquire))
        return g_global_seed.load(std::memory_order_relaxed);

    // never set: choose one from entropy, once, however many threads ask at the same time
    std::lock_guard<std::mutex> lock(g_global_seed_mutex);
    if (!g_global_seed_chosen.load(std::memory_order_relaxed))
    {
        g_global_seed.store(EntropySeed(), std::memory_order_relaxed);
        g_global_seed_chosen.store(true, std::memory_order_release);
    }
    return g_global_seed.load(std::memory_order_relaxed);
}

unsigned long SetGlobalSeedFromEntropy()
{
    unsigned long const seed = EntropySeed();
    SetGlobalSeed(seed);
    return seed;
}

unsigned long ApplyClassicRandomSeed(unsigned long randomseed)
{
    if (randomseed == 0)
        return SetGlobalSeedFromEntropy();
    SetGlobalSeed(randomseed);
    return randomseed;
}

void SetGlobalSeed(unsigned long seed)
{
    g_global_seed.store(seed, std::memory_order_relaxed);
    g_global_seed_chosen.store(true, std::memory_order_release);
    g_solve_ordinal.store(0, std::memory_order_relaxed);   // (master, ordinal) restarts here
    // the setup stream: domain = setup, index = 0.  Path/worker streams use other domains, so none
    // of them can ever reproduce this stream (the old ReseedThisThread(0) == SetGlobalSeed collision).
    // The pinned stream (b2rand/1, ADR-0044) is the source of every draw.
    records::ThreadDrawStream().Reseed(static_cast<uint64_t>(seed), kDomainSetup, 0);
}

void ReseedThisThread(uint64_t stream_key)
{
    // per-path / per-thread stream: domain = path, index = stream_key.  Deterministic from the global
    // seed and distinct for every stream_key (and distinct from the setup stream).
    records::ThreadDrawStream().Reseed(static_cast<uint64_t>(GetGlobalSeed()), kDomainPath, stream_key);
}

// Derive a distinct, deterministic child seed for a worker rank, from the master seed.  The manager
// computes these and hands one to each worker, which then calls SetGlobalSeed(child) -- so every
// process has its own non-overlapping deterministic stream, all reproducible from the one user seed.
unsigned long DerivedWorkerSeed(uint64_t worker_index)
{
    uint64_t s = static_cast<uint64_t>(GetGlobalSeed());
    uint64_t h = splitmix64(s ^ kDomainWorker ^ splitmix64(worker_index));
    if (h == 0) h = 1; // classic input reads randomseed 0 as "draw from entropy"; never hand out 0
    return static_cast<unsigned long>(h);
}

unsigned long AdvanceSolveSeed()
{
    // (master, solve ordinal) -> child seed, pure integer arithmetic (splitmix64):
    // identical on every platform, and DELIBERATELY independent of the thread-local
    // draw streams -- after a multithreaded solve the calling thread's stream position
    // depends on which paths it happened to track, so deriving from a stream would be
    // scheduling-dependent.  Since the solve surface rekeys to the returned seed
    // (SetGlobalSeed resets the ordinal), consecutive seedless solves form a
    // deterministic seed chain from the initial master.  Masked to 32 bits so the
    // value is identical on LP64 and LLP64 (Windows) platforms.
    uint64_t const master = static_cast<uint64_t>(GetGlobalSeed());
    uint64_t const ordinal = g_solve_ordinal.fetch_add(1, std::memory_order_relaxed);
    uint64_t const h = splitmix64(master ^ kDomainSolve ^ splitmix64(ordinal));
    unsigned long s = static_cast<unsigned long>(h & 0xFFFFFFFFull);
    if (s == 0) s = static_cast<unsigned long>(h >> 32);   // deterministic nonzero fallback
    if (s == 0) s = 1;                                     // so a recorded seed replays through a classic
                                                           // input file, where randomseed 0 means entropy
    return s;
}

namespace {

constexpr char RandomStateVersion[] = "b2randstate/1";

// the rest of `line` after `label` and one space, or false if the line is not labelled so
bool AfterLabel(std::string const& line, std::string const& label, std::string& rest)
{
    if (line.compare(0, label.size() + 1, label + ' ') != 0)
        return false;
    rest = line.substr(label.size() + 1);
    return true;
}

} // anon namespace

std::string RandomState::ToText() const
{
    std::ostringstream out;
    out << RandomStateVersion << '\n'
        << "global_seed " << global_seed << '\n'
        << "seed_chosen " << (seed_chosen ? 1 : 0) << '\n'
        << "solve_ordinal " << solve_ordinal << '\n'
        << "draw_stream " << draw_stream.ToText() << '\n';
    return out.str();
}

RandomState RandomState::FromText(std::string const& text)
{
    auto refuse = [](std::string const& why) {
        throw std::invalid_argument(std::string("not a ") + RandomStateVersion
                                    + " random state: " + why);
    };
    std::istringstream in(text);
    std::string line, rest;
    if (!std::getline(in, line) || line != RandomStateVersion)
        refuse("the first line is \"" + line + "\"");

    RandomState state;
    int chosen;
    if (!std::getline(in, line) || !AfterLabel(line, "global_seed", rest)
        || !(std::istringstream(rest) >> state.global_seed))
        refuse("expected the line \"global_seed <number>\", read \"" + line + "\"");
    if (!std::getline(in, line) || !AfterLabel(line, "seed_chosen", rest)
        || !(std::istringstream(rest) >> chosen) || (chosen != 0 && chosen != 1))
        refuse("expected the line \"seed_chosen <0 or 1>\", read \"" + line + "\"");
    state.seed_chosen = (chosen == 1);
    if (!std::getline(in, line) || !AfterLabel(line, "solve_ordinal", rest)
        || !(std::istringstream(rest) >> state.solve_ordinal))
        refuse("expected the line \"solve_ordinal <number>\", read \"" + line + "\"");

    if (!std::getline(in, line) || !AfterLabel(line, "draw_stream", rest))
        refuse("expected the line \"draw_stream ...\", read \"" + line + "\"");
    try
    {
        state.draw_stream = records::DrawStream::FromText(rest);
    }
    catch (std::invalid_argument const& e)
    {
        refuse(e.what());
    }

    if (std::getline(in, line) && !line.empty())
        refuse("unexpected text after the draw_stream line: \"" + line + "\"");
    return state;
}

bool RandomState::operator==(RandomState const& other) const
{
    return global_seed == other.global_seed && seed_chosen == other.seed_chosen
           && solve_ordinal == other.solve_ordinal && draw_stream == other.draw_stream;
}

RandomState GetRandomState()
{
    RandomState state;
    state.global_seed = g_global_seed.load(std::memory_order_relaxed);
    state.seed_chosen = g_global_seed_chosen.load(std::memory_order_acquire);
    state.solve_ordinal = g_solve_ordinal.load(std::memory_order_relaxed);
    state.draw_stream = records::ThreadDrawStream();
    return state;
}

void SetRandomState(RandomState const& state)
{
    g_global_seed.store(state.global_seed, std::memory_order_relaxed);
    g_global_seed_chosen.store(state.seed_chosen, std::memory_order_release);
    g_solve_ordinal.store(state.solve_ordinal, std::memory_order_relaxed);
    records::ThreadDrawStream() = state.draw_stream;
}



    real_mp RandomMp()
    {
        // ThreadPrecision (thread-local) rather than DefaultPrecision (global):
        // random numbers are generated during tracking, which may run on a
        // std::thread worker whose precision differs from the global default.
        return RandomMp(bertini::ThreadPrecision());
    }

    real_mp RandomMp(unsigned num_digits)
    {

        real_mp a;
        if (num_digits<=50)
            a = RandomMp<50>();
        else if (num_digits<=100)
            a = RandomMp<100>();
        else if (num_digits<=200)
            a = RandomMp<200>();
        else if (num_digits<=400)
            a = RandomMp<400>();
        else if (num_digits<=800)
            a = RandomMp<800>();
        else if (num_digits<=1600)
            a = RandomMp<1600>();
        else if (num_digits<=3200)
            a = RandomMp<3200>();
        else if (num_digits<=6400)
            a = RandomMp<6400>();
        else if (num_digits<=8000)
            a = RandomMp<8000>();
        else if (num_digits<=10000)
            a = RandomMp<10000>();
        else if (num_digits<=12000)
            a = RandomMp<12000>();
        else if (num_digits<=14000)
            a = RandomMp<14000>();
        else if (num_digits<=16000)
            a = RandomMp<16000>();
        else if (num_digits<=18000)
            a = RandomMp<18000>();
        else if (num_digits<=20000)
            a = RandomMp<20000>();
        else if (num_digits<=40000)
            a = RandomMp<40000>();
        else
            throw std::out_of_range("requesting random long number of digits -- higher than 40000.  this throw can be remedied by adding more cases to the generating function RandomMp in random.cpp.  If you have a better solution to this problem, please write the authors of this software.");
        a.precision(num_digits);
        return a;
    }


    void RandomMpAssign(real_mp & a, unsigned num_digits)
    {


        real_mp temp;
        temp = RandomMp(num_digits);
        a = std::move(temp);
    }




    real_mp RandomMp(const real_mp & a, const real_mp & b)
    {
        // see RandomMp() above for why ThreadPrecision rather than DefaultPrecision
        return RandomMp(a,b,bertini::ThreadPrecision());
    }

    real_mp RandomMp(const real_mp & a, const real_mp & b, unsigned num_digits)
    {

        real_mp result;
        if (num_digits<=50)
            result = RandomMp<50>(a,b);
        else if (num_digits<=100)
            result = RandomMp<100>(a,b);
        else if (num_digits<=200)
            result = RandomMp<200>(a,b);
        else if (num_digits<=400)
            result = RandomMp<400>(a,b);
        else if (num_digits<=800)
            result = RandomMp<800>(a,b);
        else if (num_digits<=1600)
            result = RandomMp<1600>(a,b);
        else if (num_digits<=3200)
            result = RandomMp<3200>(a,b);
        else if (num_digits<=6400)
            result = RandomMp<6400>(a,b);
        else if (num_digits<=8000)
            result = RandomMp<8000>(a,b);
        else if (num_digits<=10000)
            result = RandomMp<10000>(a,b);
        else if (num_digits<=12000)
            result = RandomMp<12000>(a,b);
        else if (num_digits<=14000)
            result = RandomMp<14000>(a,b);
        else if (num_digits<=16000)
            result = RandomMp<16000>(a,b);
        else if (num_digits<=18000)
            result = RandomMp<18000>(a,b);
        else if (num_digits<=20000)
            result = RandomMp<20000>(a,b);
        else if (num_digits<=40000)
            result = RandomMp<40000>(a,b);
        else
            throw std::out_of_range("requesting random long number of digits -- higher than 40000.  this throw can be remedied by adding more cases to the generating function RandomMp in random.cpp.  If you have a better solution to this problem, please write the authors of this software.");
        result.precision(num_digits);
        return result;
    }









} // namespace bertini
