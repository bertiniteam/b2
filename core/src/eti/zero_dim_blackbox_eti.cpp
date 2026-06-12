//This file is part of Bertini 2.
//
//src/eti/zero_dim_blackbox_eti.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/eti/zero_dim_blackbox_eti.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/eti/zero_dim_blackbox_eti.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file zero_dim_blackbox_eti.cpp

Explicit instantiation definitions for the ZeroDim combos reachable through
the blackbox switch ladder beyond the TotalDegree six (see zero_dim_eti.cpp):
MHomogeneous start (CloneGiven) and User start (RefToGiven), per
blackbox/config.hpp's StorageSelector.  Pairs with the extern block at the
bottom of bertini2/nag_algorithms/zero_dim_solve.hpp.  See ADR-0014.
*/

#include "bertini2/nag_algorithms/zero_dim_solve.hpp"
#include "bertini2/endgames.hpp"
#include "bertini2/system/start_systems.hpp"

namespace bertini{ namespace algorithm{

using DPT  = tracking::DoublePrecisionTracker;
using MPT  = tracking::MultiplePrecisionTracker;
using AMPT = tracking::AMPTracker;

template struct ZeroDim<DPT,  typename endgame::EndgameSelector<DPT>::PSEG,    System, start_system::MHomogeneous>;
template struct ZeroDim<DPT,  typename endgame::EndgameSelector<DPT>::Cauchy,  System, start_system::MHomogeneous>;
template struct ZeroDim<MPT,  typename endgame::EndgameSelector<MPT>::PSEG,    System, start_system::MHomogeneous>;
template struct ZeroDim<MPT,  typename endgame::EndgameSelector<MPT>::Cauchy,  System, start_system::MHomogeneous>;
template struct ZeroDim<AMPT, typename endgame::EndgameSelector<AMPT>::PSEG,   System, start_system::MHomogeneous>;
template struct ZeroDim<AMPT, typename endgame::EndgameSelector<AMPT>::Cauchy, System, start_system::MHomogeneous>;

template struct ZeroDim<DPT,  typename endgame::EndgameSelector<DPT>::PSEG,    System, start_system::User, policy::RefToGiven>;
template struct ZeroDim<DPT,  typename endgame::EndgameSelector<DPT>::Cauchy,  System, start_system::User, policy::RefToGiven>;
template struct ZeroDim<MPT,  typename endgame::EndgameSelector<MPT>::PSEG,    System, start_system::User, policy::RefToGiven>;
template struct ZeroDim<MPT,  typename endgame::EndgameSelector<MPT>::Cauchy,  System, start_system::User, policy::RefToGiven>;
template struct ZeroDim<AMPT, typename endgame::EndgameSelector<AMPT>::PSEG,   System, start_system::User, policy::RefToGiven>;
template struct ZeroDim<AMPT, typename endgame::EndgameSelector<AMPT>::Cauchy, System, start_system::User, policy::RefToGiven>;

}} // namespaces
