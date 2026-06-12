//This file is part of Bertini 2.
//
//src/eti/endgames_eti.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//src/eti/endgames_eti.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with src/eti/endgames_eti.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file endgames_eti.cpp

Explicit instantiation definitions for the closed universe of endgame types:
{PowerSeries, Cauchy} x {fixed double, fixed multiple, AMP}.  These pair with
the extern template declarations at the bottom of bertini2/endgames.hpp and
prevent every consumer TU (python bindings, blackbox, tests) from re-emitting
the entire endgame+tracker instantiation cone.  See ADR-0014.

Adding a new endgame flavor or tracker?  Add its combinations here AND to the
extern block in endgames.hpp.
*/

#include "bertini2/endgames.hpp"

namespace bertini{ namespace endgame{

using DPT = tracking::DoublePrecisionTracker;
using MPT = tracking::MultiplePrecisionTracker;

// the bases first: explicitly instantiating a derived class does not
// instantiate the base's members.
template class EndgameBase<PowerSeriesEndgame<FixedPrecEndgame<DPT>>, FixedPrecEndgame<DPT>>;
template class EndgameBase<PowerSeriesEndgame<FixedPrecEndgame<MPT>>, FixedPrecEndgame<MPT>>;
template class EndgameBase<PowerSeriesEndgame<AMPEndgame>, AMPEndgame>;
template class EndgameBase<CauchyEndgame<FixedPrecEndgame<DPT>>, FixedPrecEndgame<DPT>>;
template class EndgameBase<CauchyEndgame<FixedPrecEndgame<MPT>>, FixedPrecEndgame<MPT>>;
template class EndgameBase<CauchyEndgame<AMPEndgame>, AMPEndgame>;

template class PowerSeriesEndgame<FixedPrecEndgame<DPT>>;
template class PowerSeriesEndgame<FixedPrecEndgame<MPT>>;
template class PowerSeriesEndgame<AMPEndgame>;
template class CauchyEndgame<FixedPrecEndgame<DPT>>;
template class CauchyEndgame<FixedPrecEndgame<MPT>>;
template class CauchyEndgame<AMPEndgame>;

}} // namespaces
