#pragma once


#include "bertini2/system.hpp"
#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/endgames.hpp"

namespace bertini{
namespace blackbox{

namespace type{
enum class Start{ TotalDegree, MHom, User};
enum class Tracker{ FixedDouble, FixedMultiple, Adaptive};
enum class Endgame{ PowerSeries, Cauchy};
}

struct ZeroDimRT
{
	type::Start start = type::Start::TotalDegree;
	type::Tracker tracker = type::Tracker::Adaptive;
	type::Endgame endgame = type::Endgame::Cauchy;
};


template <typename StartType, typename TrackerType, typename EndgameType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecify(ConstTs const& ...ts)
{
	return std::make_unique<
			algorithm::ZeroDim<
				TrackerType, 
				EndgameType, 
				System, 
				StartType>
			>(ts...);

}

template <typename StartType, typename TrackerType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyEndgame(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.endgame)
	{
		case type::Endgame::PowerSeries:
			return ZeroDimSpecify<StartType, TrackerType, 
					typename tracking::EndgameSelector<TrackerType>::PSEG>(ts...);

		case type::Endgame::Cauchy:
			return ZeroDimSpecify<StartType, TrackerType, 
					typename tracking::EndgameSelector<TrackerType>::Cauchy>(ts...);
	}
}

template <typename StartType, typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyTracker(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.tracker)
	{
		case type::Tracker::FixedDouble:
			return ZeroDimSpecifyEndgame<StartType, tracking::DoublePrecisionTracker>(rt, ts...);
		case type::Tracker::FixedMultiple:
			return ZeroDimSpecifyEndgame<StartType, tracking::MultiplePrecisionTracker>(rt, ts...);
		case type::Tracker::Adaptive:
			return ZeroDimSpecifyEndgame<StartType, tracking::AMPTracker>(rt, ts...);
	}
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> ZeroDimSpecifyStart(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	switch (rt.start)
	{
		case type::Start::TotalDegree:
			return ZeroDimSpecifyTracker<start_system::TotalDegree>(rt, ts...);
		case type::Start::MHom:
			return ZeroDimSpecifyTracker<start_system::MHomogeneous>(rt, ts...);
		case type::Start::User:
			return ZeroDimSpecifyTracker<start_system::User>(rt, ts...);
	}
}

template <typename ... ConstTs>
std::unique_ptr<algorithm::AnyZeroDim> MakeZeroDim(ZeroDimRT const& rt, ConstTs const& ...ts)
{
	return ZeroDimSpecifyStart(rt, ts...);
}


} //ns blackbox
} //ns bertini
