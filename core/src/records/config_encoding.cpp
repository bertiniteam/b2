//This file is part of Bertini 2.
//
//config_encoding.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//config_encoding.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with config_encoding.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file config_encoding.cpp

\brief The canonical config encoders (ADR-0043).  Every emitted byte is part of the
persistent digest contract: changes are digest-breaking and require a `b2cfgenc`
version bump plus a golden-fixture update in the same commit.
*/

#include "bertini2/records/config_encoding.hpp"

#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <stdexcept>

#include <boost/json.hpp>

namespace bertini {
namespace records {

namespace {

	// A double as its exact IEEE-754 bit pattern: d64:<16 lowercase hex, big-endian>.
	// Bit patterns are the only encoding of a double that is single-valued everywhere;
	// decimal and %a-style hexfloat renderings vary by implementation.
	std::string ExactDouble(double v)
	{
		static_assert(sizeof(double) == sizeof(std::uint64_t), "IEEE-754 binary64 expected");
		std::uint64_t bits = 0;
		std::memcpy(&bits, &v, sizeof(bits));
		static constexpr char digits[] = "0123456789abcdef";
		std::string out = "d64:";
		for (int shift = 60; shift >= 0; shift -= 4)
			out.push_back(digits[(bits >> shift) & 0xf]);
		return out;
	}

	std::string ExactRational(mpq_rational const& q)
	{
		return q.str();
	}

	// A name as a netstring (matches the node encoder's convention, ADR-0042).
	std::string Netstring(std::string const& s)
	{
		std::ostringstream out;
		out << s.size() << ':' << s;
		return out.str();
	}

} // unnamed namespace


std::string CanonicalName(tracking::Predictor p)
{
	using tracking::Predictor;
	switch (p)
	{
		case Predictor::Constant: return "Constant";
		case Predictor::Euler: return "Euler";
		case Predictor::Heun: return "Heun";
		case Predictor::RK4: return "RK4";
		case Predictor::HeunEuler: return "HeunEuler";
		case Predictor::RKNorsett34: return "RKNorsett34";
		case Predictor::RKF45: return "RKF45";
		case Predictor::RKCashKarp45: return "RKCashKarp45";
		case Predictor::RKDormandPrince56: return "RKDormandPrince56";
		case Predictor::RKVerner67: return "RKVerner67";
	}
	throw std::runtime_error("unencodable Predictor");
}

std::string CanonicalName(algorithm::classic::AlgoChoice a)
{
	using algorithm::classic::AlgoChoice;
	switch (a)
	{
		case AlgoChoice::EvalFunctions: return "EvalFunctions";
		case AlgoChoice::EvalFunctionJacobian: return "EvalFunctionJacobian";
		case AlgoChoice::NewtonIteration: return "NewtonIteration";
		case AlgoChoice::NewtonIterationCondNum: return "NewtonIterationCondNum";
		case AlgoChoice::ZeroDim: return "ZeroDim";
		case AlgoChoice::NID: return "NID";
		case AlgoChoice::SampleComponent: return "SampleComponent";
		case AlgoChoice::MembershipTest: return "MembershipTest";
		case AlgoChoice::ExtractWitnessSet: return "ExtractWitnessSet";
		case AlgoChoice::WitnessSetProjection: return "WitnessSetProjection";
		case AlgoChoice::IsosingularStab: return "IsosingularStab";
	}
	throw std::runtime_error("unencodable AlgoChoice");
}

std::string CanonicalName(algorithm::classic::EndgameChoice e)
{
	using algorithm::classic::EndgameChoice;
	switch (e)
	{
		case EndgameChoice::PowerSeries: return "PowerSeries";
		case EndgameChoice::Cauchy: return "Cauchy";
	}
	throw std::runtime_error("unencodable EndgameChoice");
}


std::string CanonicalEncoding(tracking::SteppingConfig const& c)
{
	std::ostringstream out;
	out << "(cfg Stepping"
	    << " initial_step_size=" << ExactRational(c.initial_step_size)
	    << " max_step_size=" << ExactRational(c.max_step_size)
	    << " min_step_size=" << ExactDouble(c.min_step_size)
	    << " step_size_success_factor=" << ExactRational(c.step_size_success_factor)
	    << " step_size_fail_factor=" << ExactRational(c.step_size_fail_factor)
	    << " consecutive_successful_steps_before_stepsize_increase="
	        << c.consecutive_successful_steps_before_stepsize_increase
	    << " min_num_steps=" << c.min_num_steps
	    << " max_num_steps=" << c.max_num_steps
	    << " frequency_of_CN_estimation=" << c.frequency_of_CN_estimation
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(tracking::NewtonConfig const& c)
{
	std::ostringstream out;
	out << "(cfg Newton"
	    << " max_num_newton_iterations=" << c.max_num_newton_iterations
	    << " min_num_newton_iterations=" << c.min_num_newton_iterations
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(tracking::FixedPrecisionConfig const& c)
{
	std::ostringstream out;
	out << "(cfg FixedPrecision precision=" << c.precision << ")";
	return out.str();
}

std::string CanonicalEncoding(tracking::AdaptiveMultiplePrecisionConfig const& c)
{
	std::ostringstream out;
	out << "(cfg AdaptiveMultiplePrecision"
	    << " coefficient_bound=" << ExactDouble(c.coefficient_bound)
	    << " degree_bound=" << ExactDouble(c.degree_bound)
	    << " epsilon=" << ExactDouble(c.epsilon)
	    << " Phi=" << ExactDouble(c.Phi)
	    << " Psi=" << ExactDouble(c.Psi)
	    << " safety_digits_1=" << c.safety_digits_1
	    << " safety_digits_2=" << c.safety_digits_2
	    << " maximum_precision=" << c.maximum_precision
	    << " max_num_precision_decreases=" << c.max_num_precision_decreases
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(tracking::Predictor p)
{
	return "(cfg Predictor choice=" + CanonicalName(p) + ")";
}


std::string CanonicalEncoding(endgame::SecurityConfig const& c)
{
	std::ostringstream out;
	out << "(cfg Security"
	    << " level=" << c.level
	    << " max_norm=" << ExactDouble(c.max_norm)
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(endgame::EndgameConfig const& c)
{
	std::ostringstream out;
	out << "(cfg Endgame"
	    << " sample_point_refinement_factor=" << ExactDouble(c.sample_point_refinement_factor)
	    << " num_sample_points=" << c.num_sample_points
	    << " min_track_time=" << ExactDouble(c.min_track_time)
	    << " sample_factor=" << ExactRational(c.sample_factor)
	    << " max_num_refinements=" << c.max_num_refinements
	    << " final_tolerance=" << ExactDouble(c.final_tolerance)
	    << " refine_when_increasing_precision=" << (c.refine_when_increasing_precision ? 1 : 0)
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(endgame::PowerSeriesConfig const& c)
{
	std::ostringstream out;
	out << "(cfg PowerSeries"
	    << " max_cycle_number=" << c.max_cycle_number
	    << " cycle_number_amplification=" << c.cycle_number_amplification
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(endgame::CauchyConfig const& c)
{
	std::ostringstream out;
	out << "(cfg Cauchy"
	    << " cycle_cutoff_time=" << ExactDouble(c.cycle_cutoff_time)
	    << " ratio_cutoff_time=" << ExactDouble(c.ratio_cutoff_time)
	    << " minimum_for_c_over_k_stabilization=" << ExactDouble(c.minimum_for_c_over_k_stabilization)
	    << " num_needed_for_stabilization=" << c.num_needed_for_stabilization
	    << " maximum_cauchy_ratio=" << ExactDouble(c.maximum_cauchy_ratio)
	    << " fail_safe_maximum_cycle_number=" << c.fail_safe_maximum_cycle_number
	    << " num_consecutive_same_cycle_number=" << c.num_consecutive_same_cycle_number
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(endgame::TrackBackConfig const& c)
{
	std::ostringstream out;
	out << "(cfg TrackBack"
	    << " minimum_cycle=" << c.minimum_cycle
	    << " junk_removal_test=" << (c.junk_removal_test ? 1 : 0)
	    << " max_depth_LDT=" << c.max_depth_LDT
	    << ")";
	return out.str();
}


std::string CanonicalEncoding(algorithm::TolerancesConfig const& c)
{
	std::ostringstream out;
	out << "(cfg Tolerances"
	    << " newton_before_endgame=" << ExactDouble(c.newton_before_endgame)
	    << " newton_during_endgame=" << ExactDouble(c.newton_during_endgame)
	    << " final_tolerance=" << ExactDouble(c.final_tolerance)
	    << " path_truncation_threshold=" << ExactDouble(c.path_truncation_threshold)
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(algorithm::MidPathConfig const& c)
{
	std::ostringstream out;
	out << "(cfg MidPath same_point_tolerance=" << ExactDouble(c.same_point_tolerance) << ")";
	return out.str();
}

std::string CanonicalEncoding(algorithm::AutoRetrackConfig const& c)
{
	std::ostringstream out;
	out << "(cfg AutoRetrack midpath_decrease_tolerance_factor="
	    << ExactDouble(c.midpath_decrease_tolerance_factor) << ")";
	return out.str();
}

std::string CanonicalEncoding(algorithm::SharpeningConfig const& c)
{
	std::ostringstream out;
	out << "(cfg Sharpening"
	    << " sharpendigits=" << c.sharpendigits
	    << " function_residual_tolerance=" << ExactDouble(c.function_residual_tolerance)
	    << " ratio_tolerance=" << ExactDouble(c.ratio_tolerance)
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(algorithm::RegenerationConfig const& c)
{
	std::ostringstream out;
	out << "(cfg Regeneration"
	    << " remove_infinite_endpoints=" << (c.remove_infinite_endpoints ? 1 : 0)
	    << " higher_dimension_check=" << (c.higher_dimension_check ? 1 : 0)
	    << " start_level=" << c.start_level
	    << " slice_newton_before_endgame=" << ExactDouble(c.slice_newton_before_endgame)
	    << " slice_newton_during_endgame=" << ExactDouble(c.slice_newton_during_endgame)
	    << " slice_final_tolerance=" << ExactDouble(c.slice_final_tolerance)
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(algorithm::PostProcessingConfig const& c)
{
	std::ostringstream out;
	out << "(cfg PostProcessing"
	    << " real_threshold=" << ExactDouble(c.real_threshold)
	    << " endpoint_finite_threshold=" << ExactDouble(c.endpoint_finite_threshold)
	    << " same_point_tolerance_multiplier=" << ExactDouble(c.same_point_tolerance_multiplier)
	    << " condition_number_threshold=" << ExactDouble(c.condition_number_threshold)
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(algorithm::ZeroDimConfig const& c)
{
	// num_threads deliberately excluded: transient (must not change what was computed).
	std::ostringstream out;
	out << "(cfg ZeroDim"
	    << " initial_ambient_precision=" << c.initial_ambient_precision
	    << " max_num_crossed_path_resolve_attempts=" << c.max_num_crossed_path_resolve_attempts
	    << " start_time=" << ExactRational(c.start_time)
	    << " endgame_boundary=" << ExactRational(c.endgame_boundary)
	    << " target_time=" << ExactRational(c.target_time)
	    << " path_variable_name=" << Netstring(c.path_variable_name)
	    << ")";
	return out.str();
}

std::string CanonicalEncoding(algorithm::MetaConfig const& c)
{
	return "(cfg Meta tracktype=" + CanonicalName(c.tracktype) + ")";
}

std::string CanonicalEncoding(algorithm::classic::EndgameChoiceConfig const& c)
{
	return "(cfg EndgameChoice endgame=" + CanonicalName(c.endgame) + ")";
}

namespace {

	// The canonical text's exact scalar encodings, decoded for the derived JSON view:
	// d64 bit patterns become shortest round-trip decimals (a reader sees 1e-05, not
	// "d64:3ee4f8b588e368f1"; the shortest form parses back to the identical double),
	// netstrings shed their length prefix, and plain integers become JSON numbers.
	// The DIGEST contract stays on the canonical text; this only affects presentation.
	std::string ValueAsJsonToken(std::string const& v)
	{
		if (v.rfind("d64:", 0) == 0 && v.size() == 20)
		{
			std::uint64_t bits = 0;
			bool ok = true;
			for (std::size_t i = 4; i < v.size(); ++i)
			{
				auto const c = v[i];
				bits <<= 4;
				if (c >= '0' && c <= '9') bits |= static_cast<std::uint64_t>(c - '0');
				else if (c >= 'a' && c <= 'f') bits |= static_cast<std::uint64_t>(c - 'a' + 10);
				else { ok = false; break; }
			}
			if (ok)
			{
				double d = 0;
				std::memcpy(&d, &bits, sizeof(d));
				if (std::isfinite(d))
				{
					char buffer[32];
					auto const res = std::to_chars(buffer, buffer + sizeof(buffer), d);
					return std::string(buffer, res.ptr);   // shortest round-trip decimal
				}
			}
		}
		if (auto const colon = v.find(':'); colon != std::string::npos)
		{
			// netstring "N:payload" -> the payload, JSON-escaped
			bool numeric_prefix = colon > 0;
			for (std::size_t i = 0; i < colon && numeric_prefix; ++i)
				numeric_prefix = (v[i] >= '0' && v[i] <= '9');
			if (numeric_prefix
			    && std::stoull(v.substr(0, colon)) == v.size() - colon - 1)
				return boost::json::serialize(boost::json::value(v.substr(colon + 1)));
		}
		if (!v.empty() && v.find_first_not_of("-0123456789") == std::string::npos
		    && v.find('-', 1) == std::string::npos && v.size() <= 18)
			return v;   // a plain integer: emit as a JSON number
		return boost::json::serialize(boost::json::value(v));   // exact string (e.g. "1/10")
	}

} // unnamed namespace

std::string ConfigTextAsJson(std::string const& canonical_text, std::string const& digest_hex)
{
	std::ostringstream out;
	out << "{\n \"schema\": \"" << ConfigEncodingVersion << "\",\n"
	    << " \"digest\": \"" << digest_hex << "\",\n"
	    << " \"configs\": {";
	bool first_config = true;
	std::istringstream lines(canonical_text);
	for (std::string line; std::getline(lines, line); )
	{
		if (line.rfind("(cfg ", 0) != 0)
			continue;
		// "(cfg Name k=v k=v ...)": split the name, then the k=v fields
		std::string const body = line.substr(5, line.size() - 6);   // drop "(cfg " and ")"
		auto const name_end = body.find(' ');
		std::string const name = body.substr(0, name_end);
		out << (first_config ? "" : ",") << "\n  \"" << name << "\": {";
		first_config = false;
		bool first_field = true;
		if (name_end != std::string::npos)
		{
			std::istringstream fields(body.substr(name_end + 1));
			for (std::string field; std::getline(fields, field, ' '); )
			{
				auto const eq = field.find('=');
				if (eq == std::string::npos)
					continue;
				out << (first_field ? "" : ",") << "\n   \"" << field.substr(0, eq)
				    << "\": " << ValueAsJsonToken(field.substr(eq + 1));
				first_field = false;
			}
		}
		out << "\n  }";
	}
	out << "\n }\n}\n";
	return out.str();
}

} // namespace records
} // namespace bertini
