//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/output.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/output.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/output.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/nag_algorithms/output.hpp 

\brief Provides some outputting classes for printing results of running algorithms

You can also provide your own, that's the point of these policies.
*/

#pragma once

#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/io/generators.hpp"


namespace bertini {
namespace algorithm {
namespace output {

/**
\brief Prints data from algorithms
*/
template <typename AlgoT>
struct Outputter
{ };


template <typename ...T>
struct Classic
{};


template <typename A, typename B, typename C, typename D, template<typename, typename> class E>
struct Classic <ZeroDim<A,B,C,D,E>>
{
	using ZDT = ZeroDim<A,B,C,D,E>;

	template <typename OutT>
	static
	void All(OutT & out, ZDT const& zd)
	{
		MainData(out, zd);
	}

	
	template <typename OutT>
	static
	void MainData(OutT & out, ZDT const& zd)
	{
		Variables(out, zd);
		out << '\n';

		const auto& s = zd.FinalSolutions();
		const auto& m = zd.FinalSolutionMetadata();
		const auto n = s.size();
		for (decltype(s.size()) ii=0; ii<n; ++ii)
		{
			EndPointMD(ii, out, zd);
			out << '\n';
			EndPoint(ii, out, zd);
			out << '\n';
		}
	}

	template <typename OutT>
	static 
	void Variables(OutT & out, ZDT const& zd)
	{
		const auto& sys = zd.TargetSystem();

		out << sys.NumVariables() << "\n\n";
		const auto& vars = sys.Variables();
		for (const auto& x : vars)
			out << *x << ' ';
		out << '\n';
	}

	template <typename IndexT, typename OutT>
	static
	void EndPoint(IndexT const& ind, OutT & out, ZDT const& zd)
	{
		generators::Classic::generate(boost::spirit::ostream_iterator(out), zd.FinalSolutions()[ind]);
	}

	template <typename IndexT, typename OutT>
	static
	void EndPointMD(IndexT const& ind, OutT & out, ZDT const& zd)
	{
		const auto& data = zd.FinalSolutionMetadata()[ind];
		out << data.path_index << '\n'
			<< data.solution_index << '\n'
			<< data.function_residual << '\n'
			<< data.condition_number << '\n'
			<< data.newton_residual << '\n';
		generators::Classic::generate(boost::spirit::ostream_iterator(out), data.final_time_used);
			// << data.final_time_used << '\n'
		out << '\n' << data.max_precision_used << '\n'
			<< data.accuracy_estimate << '\n'
			<< data.cycle_num << '\n'
			<< data.multiplicity << '\n';	
	}
};

} // namespace output
} // namespace algorithm
} // namespace bertini
