//This file is part of Bertini 2.
//
//solver_recording.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//solver_recording.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with solver_recording.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file solver_recording.hpp

\brief FullPathResult <-> track-record serialization (ledgerrec/1; the arc's rung 4).

The whole-path result (parallel/path_result.hpp) is the manager-side unit of completed
work in every topology (serial / threaded / MPI): StoreFullPathResult is the single
installation point.  So a track record IS a serialized FullPathResult, and hydration is
decoding one and replaying it through the same installer -- hydrated solver state is
identical to computed state BY CONSTRUCTION.

Scalar encodings are exact and human-readable: doubles as %.17g decimal (round-trips
bit-exactly per IEEE-754), multiprecision values as full-precision decimal strings plus
their precision.  Points are the treasure -- they stay readable.
*/

#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>

#include <boost/json.hpp>

#include "bertini2/num_traits.hpp"
#include "bertini2/parallel/path_result.hpp"

namespace bertini {
namespace records {

/// \brief A double as full-round-trip decimal text (%.17g: bit-exact per IEEE-754).
inline std::string ExactDoubleText(double v)
{
	char buffer[32];
	std::snprintf(buffer, sizeof(buffer), "%.17g", v);
	return buffer;
}

/// \brief Parse a double written by ExactDoubleText (bit-exact round trip).
inline double DoubleFromText(std::string const& s)
{
	return std::strtod(s.c_str(), nullptr);
}

/// \brief Encode one real scalar as text: doubles via %.17g; real_mp at full stored precision.
inline std::string EncodeRealScalar(double v) { return ExactDoubleText(v); }
/// \brief Encode one multiprecision real scalar as full-precision decimal text.
inline std::string EncodeRealScalar(real_mp const& v) { return v.str(0, std::ios::scientific); }

/// \brief Encode one complex value as ["re", "im", precision] (precision in digits).
inline boost::json::array EncodeComplexScalar(complex_dbl const& z)
{
	return {ExactDoubleText(z.real()), ExactDoubleText(z.imag()),
	        static_cast<std::int64_t>(DoublePrecision())};
}
/// \brief Encode one multiprecision complex value as ["re", "im", precision].
inline boost::json::array EncodeComplexScalar(complex_mp const& z)
{
	return {z.real().str(0, std::ios::scientific), z.imag().str(0, std::ios::scientific),
	        static_cast<std::int64_t>(z.precision())};
}

/// \brief Decode a complex value written by EncodeComplexScalar, at its recorded precision.
template <typename ComplexT>
ComplexT DecodeComplexScalar(boost::json::array const& triple);

/// \brief Decode a double-precision complex value.
template <>
inline complex_dbl DecodeComplexScalar<complex_dbl>(boost::json::array const& triple)
{
	return complex_dbl(DoubleFromText(std::string(triple.at(0).as_string())),
	                   DoubleFromText(std::string(triple.at(1).as_string())));
}

/// \brief Decode a multiprecision complex value at its recorded precision.
template <>
inline complex_mp DecodeComplexScalar<complex_mp>(boost::json::array const& triple)
{
	auto const digits = static_cast<unsigned>(triple.at(2).as_int64());
	complex_mp z(real_mp(std::string(triple.at(0).as_string()), digits),
	             real_mp(std::string(triple.at(1).as_string()), digits));
	z.precision(digits);
	return z;
}

/// \brief Encode a point (vector of complex coordinates) as an array of scalar triples.
template <typename ComplexT>
boost::json::array EncodePoint(Vec<ComplexT> const& point)
{
	boost::json::array out;
	for (Eigen::Index ii = 0; ii < point.size(); ++ii)
		out.push_back(EncodeComplexScalar(point(ii)));
	return out;
}

/// \brief Decode a point written by EncodePoint.
template <typename ComplexT>
Vec<ComplexT> DecodePoint(boost::json::array const& coords)
{
	Vec<ComplexT> out(static_cast<Eigen::Index>(coords.size()));
	for (std::size_t ii = 0; ii < coords.size(); ++ii)
		out(static_cast<Eigen::Index>(ii)) = DecodeComplexScalar<ComplexT>(coords[ii].as_array());
	return out;
}

/**
\brief Serialize one whole-path result as a ledgerrec/1 track record body.

Field names mirror parallel::FullPathResult.  SuccessCodes are recorded as integers
(their values are part of the ledgerrec contract).  The caller adds "kind"/"run"/
"index"/"status"/"start".
*/
template <typename ComplexT>
boost::json::object EncodeFullPathResult(parallel::FullPathResult<ComplexT> const& r)
{
	boost::json::object out;
	out["pre_endgame_success_code"] = static_cast<std::int64_t>(r.pre_endgame_success_code);
	out["boundary_point"] = EncodePoint(r.boundary_point);
	out["boundary_stepsize"] = EncodeRealScalar(r.boundary_stepsize);
	out["boundary_precision"] = static_cast<std::int64_t>(r.boundary_precision);
	out["endgame_success_code"] = static_cast<std::int64_t>(r.endgame_success_code);
	out["endpoint"] = EncodePoint(r.solution);
	out["function_residual"] = ExactDoubleText(r.function_residual);
	out["condition_number"] = ExactDoubleText(r.condition_number);
	out["newton_residual"] = ExactDoubleText(r.newton_residual);
	out["final_time_used"] = EncodeComplexScalar(r.final_time_used);
	out["accuracy_estimate"] = ExactDoubleText(r.accuracy_estimate);
	out["accuracy_estimate_user_coords"] = ExactDoubleText(r.accuracy_estimate_user_coords);
	out["cycle_num"] = static_cast<std::int64_t>(r.cycle_num);
	out["precision_digits"] = static_cast<std::int64_t>(r.precision_digits);
	out["accuracy_digits"] = static_cast<std::int64_t>(r.accuracy_digits);
	out["precision_changed"] = r.precision_changed;
	out["time_of_first_prec_increase"] = EncodeComplexScalar(r.time_of_first_prec_increase);
	out["max_precision_used"] = static_cast<std::int64_t>(r.max_precision_used);
	out["path_time_seconds"] = ExactDoubleText(r.path_time_seconds);
	return out;
}

/// \brief Reconstruct a whole-path result from a track record (the hydration half).
template <typename ComplexT>
parallel::FullPathResult<ComplexT> DecodeFullPathResult(boost::json::object const& rec,
                                                        std::size_t path_index)
{
	using RealT = typename parallel::FullPathResult<ComplexT>::RealT;
	parallel::FullPathResult<ComplexT> r;
	r.path_index = path_index;
	r.pre_endgame_success_code =
		static_cast<SuccessCode>(rec.at("pre_endgame_success_code").as_int64());
	r.boundary_point = DecodePoint<ComplexT>(rec.at("boundary_point").as_array());
	if constexpr (std::is_same<RealT, double>::value)
		r.boundary_stepsize = DoubleFromText(std::string(rec.at("boundary_stepsize").as_string()));
	else
		r.boundary_stepsize = RealT(std::string(rec.at("boundary_stepsize").as_string()));
	r.boundary_precision = static_cast<unsigned>(rec.at("boundary_precision").as_int64());
	r.endgame_success_code = static_cast<SuccessCode>(rec.at("endgame_success_code").as_int64());
	r.solution = DecodePoint<ComplexT>(rec.at("endpoint").as_array());
	r.function_residual = DoubleFromText(std::string(rec.at("function_residual").as_string()));
	r.condition_number = DoubleFromText(std::string(rec.at("condition_number").as_string()));
	r.newton_residual = DoubleFromText(std::string(rec.at("newton_residual").as_string()));
	r.final_time_used = DecodeComplexScalar<ComplexT>(rec.at("final_time_used").as_array());
	r.accuracy_estimate = DoubleFromText(std::string(rec.at("accuracy_estimate").as_string()));
	r.accuracy_estimate_user_coords =
		DoubleFromText(std::string(rec.at("accuracy_estimate_user_coords").as_string()));
	r.cycle_num = static_cast<unsigned>(rec.at("cycle_num").as_int64());
	r.precision_digits = static_cast<unsigned>(rec.at("precision_digits").as_int64());
	r.accuracy_digits = static_cast<unsigned>(rec.at("accuracy_digits").as_int64());
	r.precision_changed = rec.at("precision_changed").as_bool();
	r.time_of_first_prec_increase =
		DecodeComplexScalar<ComplexT>(rec.at("time_of_first_prec_increase").as_array());
	r.max_precision_used = static_cast<unsigned>(rec.at("max_precision_used").as_int64());
	r.path_time_seconds = DoubleFromText(std::string(rec.at("path_time_seconds").as_string()));
	return r;
}

} // namespace records
} // namespace bertini
