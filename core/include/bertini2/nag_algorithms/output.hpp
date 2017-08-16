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
		out << "\n\n\n MAINDATA \n\n\n";
		MainData(out, zd);

		out << "\n\n\n RAWDATA \n\n\n";
		RawData(out, zd);
	}

	
	template <typename OutT>
	static
	void MainData(OutT & out, ZDT const& zd)
	{
		NumVariables(out, zd);
		Variables(out, zd,"\n\n");

		const auto& s = zd.FinalSolutions();
		const auto& m = zd.FinalSolutionMetadata();
		const auto n = s.size();
		for (decltype(s.size()) ii{0}; ii<n; ++ii)
		{
			if (zd.FinalSolutionMetadata()[ii].endgame_success == SuccessCode::NeverStarted)
				continue;
			
			EndPointMDFull(ii, out, zd);
			EndPoint(ii, out, zd,"\n\n");
		}

		out << "\n\n";
		TargetSystem(out,zd,"\n\n");

		StartSystem(out,zd,"\n\n");

		Homotopy(out,zd,"\n\n");
	}

	template <typename OutT>
	static 
	void RawData(OutT & out, ZDT const& zd)
	{
		const auto n = zd.FinalSolutions().size();
		NumVariables(out, zd,"\n\n");
		for (decltype(zd.FinalSolutions().size()) ii{0}; ii<n; ++ii)
		{
			if (zd.FinalSolutionMetadata()[ii].endgame_success == SuccessCode::NeverStarted)
				continue;
			EndPointMDRaw(ii,out,zd,"\n\n");
		}

		out << "\n\n";
		TargetSystem(out,zd,"\n\n");
	}


	template <typename OutT>
	static 
	void NumVariables(OutT & out, ZDT const& zd, std::string const& additional = "\n")
	{
		const auto& sys = zd.TargetSystem();
		out << sys.NumVariables() << additional;
	}


	template <typename OutT>
	static 
	void Variables(OutT & out, ZDT const& zd, std::string const& additional = "\n")
	{
		const auto& sys = zd.TargetSystem();
		const auto& vars = sys.VariableOrdering();
		for (const auto& x : vars)
			out << *x << ' ';
		out << additional;
	}


	template <typename OutT>
	static 
	void TargetSystem(OutT & out, ZDT const& zd, std::string const& additional = "\n")
	{
		out << zd.TargetSystem() << additional;
	}

	template <typename OutT>
	static 
	void StartSystem(OutT & out, ZDT const& zd, std::string const& additional = "\n")
	{
		out << zd.StartSystem() << additional;
	}

	template <typename OutT>
	static 
	void Homotopy(OutT & out, ZDT const& zd, std::string const& additional = "\n")
	{
		out << zd.Homotopy() << additional;
	}

	template <typename IndexT, typename OutT>
	static
	void EndPoint(IndexT const& ind, OutT & out, ZDT const& zd, std::string const& additional = "")
	{	
		generators::Classic::generate(boost::spirit::ostream_iterator(out), zd.FinalSolutions()[ind]);
		out << additional;
	}

	template <typename IndexT, typename OutT>
	static
	void EndPointDehom(IndexT const& ind, OutT & out, ZDT const& zd, std::string const& additional = "")
	{	
		DefaultPrecision(Precision(zd.FinalSolutions()[ind]));

		generators::Classic::generate(boost::spirit::ostream_iterator(out), zd.TargetSystem().DehomogenizePoint(zd.FinalSolutions()[ind]));
		out << additional;
	}

	template <typename IndexT, typename OutT>
	static
	void EndPointMDFull(IndexT const& ind, OutT & out, ZDT const& zd, std::string const& additional = "\n")
	{
		const auto& data = zd.FinalSolutionMetadata()[ind];
		out << data.path_index << '\n'
			<< data.solution_index << '\n'
			<< data.condition_number << '\n'
			<< data.function_residual << '\n'
			<< data.newton_residual << '\n';
		generators::Classic::generate(boost::spirit::ostream_iterator(out), data.final_time_used);
		out << '\n' << data.max_precision_used << '\n';

		generators::Classic::generate(boost::spirit::ostream_iterator(out), data.time_of_first_prec_increase);
		out << '\n' << data.accuracy_estimate << '\n'
			<< data.accuracy_estimate_user_coords << '\n'
			<< data.cycle_num << '\n'
			<< data.multiplicity << '\n'
			<< data.pre_endgame_success << ' ' << data.endgame_success << '\n';
		out << additional;
	}



	template <typename IndexT, typename OutT>
	static
	void EndPointMDRaw(IndexT const& ind, OutT & out, ZDT const& zd, std::string const& additional = "\n")
	{
		const auto& pt = zd.FinalSolutions()[ind];
		const auto& data = zd.FinalSolutionMetadata()[ind];
		out << data.path_index << '\n'
			<< Precision(pt) << '\n';
		EndPoint(ind, out, zd);
		out << data.function_residual << '\n'
			<< data.condition_number << '\n'
			<< data.newton_residual << '\n';
		generators::Classic::generate(boost::spirit::ostream_iterator(out), data.final_time_used);
		out << '\n' << data.accuracy_estimate << '\n';
		generators::Classic::generate(boost::spirit::ostream_iterator(out), data.time_of_first_prec_increase);
		out << '\n' << data.cycle_num << '\n'
			<< data.endgame_success << '\n'
			<< additional;
	}


};

} // namespace output
} // namespace algorithm
} // namespace bertini
