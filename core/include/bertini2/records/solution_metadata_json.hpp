//This file is part of Bertini 2.
//
//solution_metadata_json.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//solution_metadata_json.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with solution_metadata_json.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file solution_metadata_json.hpp

\brief A solution's metadata as exact JSON, and back (bertiniteam/b2#521).

Every field of algorithm::SolutionMetaData is written, with the scalar encodings of the track
records (solver_recording.hpp): a double as its exact decimal, a multiprecision value as
full-precision decimal with its precision, a success code by its fixed name.  Reading the text
back gives metadata equal to the written one in every field, path_time_seconds included, so a
program that keeps solutions between runs keeps what b2 knew about them.  The text names its
format, `b2meta/1`, and its numeric type, and a reader refuses any other.
*/

#pragma once

#include <stdexcept>
#include <string>

#include <boost/json.hpp>

#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/records/solver_recording.hpp"

namespace bertini {
namespace records {

/// \brief The format a solution's metadata is written in; a new encoding is a new version.
constexpr char SolutionMetaDataFormat[] = "b2meta/1";

/**
\brief The success code a CanonicalName names.
\param name The name CanonicalName gave.
\return The code.
\throws std::invalid_argument if no code has that name.
*/
inline SuccessCode SuccessCodeFromName(std::string const& name)
{
    for (auto code : {SuccessCode::NeverStarted, SuccessCode::Success,
                      SuccessCode::HigherPrecisionNecessary, SuccessCode::ReduceStepSize,
                      SuccessCode::GoingToInfinity, SuccessCode::FailedToConverge,
                      SuccessCode::MatrixSolveFailure,
                      SuccessCode::MatrixSolveFailureFirstPartOfPrediction,
                      SuccessCode::MaxNumStepsTaken, SuccessCode::MaxPrecisionReached,
                      SuccessCode::MinStepSizeReached, SuccessCode::Failure,
                      SuccessCode::SingularStartPoint, SuccessCode::ExternallyTerminated,
                      SuccessCode::MinTrackTimeReached, SuccessCode::SecurityMaxNormReached,
                      SuccessCode::CycleNumTooHigh, SuccessCode::FailedToSelectPrecisionAndStepsize,
                      SuccessCode::WallClockLimitReached})
        if (CanonicalName(code) == name)
            return code;
    throw std::invalid_argument("no success code is named \"" + name + "\"");
}

/// \brief The name of a solution's numeric type in its metadata's text.
/// \return "double" or "multiprecision".
template <typename ComplexT>
inline std::string NumericTypeName();

/// \brief The double-precision name.
template <>
inline std::string NumericTypeName<complex_dbl>() { return "double"; }

/// \brief The multiprecision name.
template <>
inline std::string NumericTypeName<complex_mp>() { return "multiprecision"; }

/**
\brief A solution's metadata as a JSON object: every field, exactly.
\param m The metadata.
\return The object, beginning with its format and numeric type.
*/
template <typename ComplexT>
boost::json::object EncodeSolutionMetaData(algorithm::SolutionMetaData<ComplexT> const& m)
{
    boost::json::object out;
    out["format"] = SolutionMetaDataFormat;
    out["numeric_type"] = NumericTypeName<ComplexT>();
    out["path_index"] = static_cast<std::uint64_t>(m.path_index);
    out["solution_index"] = static_cast<std::uint64_t>(m.solution_index);
    out["precision_changed"] = m.precision_changed;
    out["time_of_first_prec_increase"] = EncodeComplexScalar(m.time_of_first_prec_increase);
    out["max_precision_used"] = static_cast<std::uint64_t>(m.max_precision_used);
    out["path_time_seconds"] = ExactDoubleText(m.path_time_seconds);
    out["num_successful_steps"] = static_cast<std::uint64_t>(m.num_successful_steps);
    out["num_failed_steps"] = static_cast<std::uint64_t>(m.num_failed_steps);
    out["wall_clock_limit_seconds"] = ExactDoubleText(m.wall_clock_limit_seconds);
    out["latest_path_point"] = EncodePoint(m.latest_path_point);
    out["pre_endgame_success_code"] = CanonicalName(m.pre_endgame_success_code);
    out["condition_number"] = ExactDoubleText(m.condition_number);
    {
        boost::json::array values;
        for (Eigen::Index ii = 0; ii < m.singular_values.size(); ++ii)
            values.push_back(boost::json::value(ExactDoubleText(m.singular_values(ii))));
        out["singular_values"] = std::move(values);
    }
    out["newton_residual"] = ExactDoubleText(m.newton_residual);
    out["final_time_used"] = EncodeComplexScalar(m.final_time_used);
    out["accuracy_estimate_internal_coords"] = ExactDoubleText(m.accuracy_estimate_internal_coords);
    out["accuracy_estimate_user_coords"] = ExactDoubleText(m.accuracy_estimate_user_coords);
    out["cycle_num"] = static_cast<std::uint64_t>(m.cycle_num);
    out["precision_digits"] = static_cast<std::uint64_t>(m.precision_digits);
    out["accuracy_digits"] = static_cast<std::uint64_t>(m.accuracy_digits);
    out["endgame_success_code"] = CanonicalName(m.endgame_success_code);
    out["function_residual"] = ExactDoubleText(m.function_residual);
    out["multiplicity"] = static_cast<std::int64_t>(m.multiplicity);
    out["multiplicity_representative"] = m.multiplicity_representative;
    out["is_real"] = m.is_real;
    out["is_finite"] = m.is_finite;
    out["is_singular"] = m.is_singular;
    out["crossing_unresolved"] = m.crossing_unresolved;
    out["is_nonsolution"] = m.is_nonsolution;
    return out;
}

/**
\brief The metadata an EncodeSolutionMetaData object describes.
\param in The object.
\return The metadata, equal to the encoded one in every field.
\throws std::invalid_argument if the object is of another format or numeric type, or lacks a
field or holds one of the wrong kind.
*/
template <typename ComplexT>
algorithm::SolutionMetaData<ComplexT> DecodeSolutionMetaData(boost::json::object const& in)
{
    auto const* format = in.if_contains("format");
    if (!format || !format->is_string() || format->as_string() != SolutionMetaDataFormat)
        throw std::invalid_argument(std::string("not a ") + SolutionMetaDataFormat
                                    + " solution metadata object");
    auto const* numeric = in.if_contains("numeric_type");
    if (!numeric || !numeric->is_string() || numeric->as_string() != NumericTypeName<ComplexT>())
        throw std::invalid_argument(std::string("solution metadata of another numeric type; this "
                                                "reader reads ") + NumericTypeName<ComplexT>());

    auto field = [&in](char const* name) -> boost::json::value const& {
        auto const* v = in.if_contains(name);
        if (!v)
            throw std::invalid_argument(std::string(SolutionMetaDataFormat)
                                        + " solution metadata lacks the field " + name);
        return *v;
    };
    auto number = [&field](char const* name) {
        return field(name).template to_number<std::uint64_t>();
    };
    auto text = [&field](char const* name) { return std::string(field(name).as_string()); };
    auto real = [&text](char const* name) { return DoubleFromText(text(name)); };

    algorithm::SolutionMetaData<ComplexT> m;
    try
    {
        m.path_index = number("path_index");
        m.solution_index = number("solution_index");
        m.precision_changed = field("precision_changed").as_bool();
        m.time_of_first_prec_increase =
            DecodeComplexScalar<ComplexT>(field("time_of_first_prec_increase").as_array());
        m.max_precision_used = static_cast<decltype(m.max_precision_used)>(number("max_precision_used"));
        m.path_time_seconds = real("path_time_seconds");
        m.num_successful_steps = static_cast<unsigned>(number("num_successful_steps"));
        m.num_failed_steps = static_cast<unsigned>(number("num_failed_steps"));
        m.wall_clock_limit_seconds = real("wall_clock_limit_seconds");
        m.latest_path_point = DecodePoint<ComplexT>(field("latest_path_point").as_array());
        m.pre_endgame_success_code = SuccessCodeFromName(text("pre_endgame_success_code"));
        m.condition_number = real("condition_number");
        auto const& values = field("singular_values").as_array();
        m.singular_values.resize(static_cast<Eigen::Index>(values.size()));
        for (std::size_t ii = 0; ii < values.size(); ++ii)
            m.singular_values(static_cast<Eigen::Index>(ii)) =
                DoubleFromText(std::string(values[ii].as_string()));
        m.newton_residual = real("newton_residual");
        m.final_time_used = DecodeComplexScalar<ComplexT>(field("final_time_used").as_array());
        m.accuracy_estimate_internal_coords = real("accuracy_estimate_internal_coords");
        m.accuracy_estimate_user_coords = real("accuracy_estimate_user_coords");
        m.cycle_num = static_cast<unsigned>(number("cycle_num"));
        m.precision_digits = static_cast<unsigned>(number("precision_digits"));
        m.accuracy_digits = static_cast<unsigned>(number("accuracy_digits"));
        m.endgame_success_code = SuccessCodeFromName(text("endgame_success_code"));
        m.function_residual = real("function_residual");
        m.multiplicity = static_cast<int>(field("multiplicity").template to_number<std::int64_t>());
        m.multiplicity_representative = field("multiplicity_representative").as_bool();
        m.is_real = field("is_real").as_bool();
        m.is_finite = field("is_finite").as_bool();
        m.is_singular = field("is_singular").as_bool();
        m.crossing_unresolved = field("crossing_unresolved").as_bool();
        m.is_nonsolution = field("is_nonsolution").as_bool();
    }
    catch (std::invalid_argument const&)
    {
        throw;
    }
    catch (std::exception const& e)
    {
        // boost::json's accessors throw system_error for a value of the wrong kind
        throw std::invalid_argument(std::string(SolutionMetaDataFormat)
                                    + " solution metadata has a field of the wrong kind: " + e.what());
    }
    return m;
}

/**
\brief A solution's metadata as JSON text.
\param m The metadata.
\return The text of EncodeSolutionMetaData's object.
*/
template <typename ComplexT>
std::string SolutionMetaDataToJson(algorithm::SolutionMetaData<ComplexT> const& m)
{
    return boost::json::serialize(EncodeSolutionMetaData(m));
}

/**
\brief The metadata a SolutionMetaDataToJson text describes.
\param text The text.
\return The metadata.
\throws std::invalid_argument if the text is not JSON, or is not such metadata (see
DecodeSolutionMetaData).
*/
template <typename ComplexT>
algorithm::SolutionMetaData<ComplexT> SolutionMetaDataFromJson(std::string const& text)
{
    boost::system::error_code error;
    auto const value = boost::json::parse(text, error);
    if (error || !value.is_object())
        throw std::invalid_argument(std::string("not a ") + SolutionMetaDataFormat
                                    + " solution metadata text: not a JSON object");
    return DecodeSolutionMetaData<ComplexT>(value.as_object());
}

} // namespace records
} // namespace bertini
