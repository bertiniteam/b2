//This file is part of Bertini 2.
//
//bertini2/endgames.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/endgames.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/endgames.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of notre dame

/**
\file bertini2/endgames.hpp 

\brief Collects the various header files which define the Bertini2 endgames.
*/


#pragma once

#include "bertini2/endgames/amp_endgame.hpp"
#include "bertini2/endgames/fixed_prec_endgame.hpp"

#include "bertini2/endgames/powerseries.hpp"
#include "bertini2/endgames/cauchy.hpp"

#include "bertini2/endgames/observers.hpp"


// Explicit instantiation declarations — suppress re-instantiation of the
// endgame stack in every including TU.  The definitions live in
// core/src/eti/endgames_eti.cpp; see ADR-0014.  Adding a flavor or tracker?
// Extend both lists.
namespace bertini{ namespace endgame{

extern template class EndgameBase<PowerSeriesEndgame<FixedPrecEndgame<tracking::DoublePrecisionTracker>>, FixedPrecEndgame<tracking::DoublePrecisionTracker>>;
extern template class EndgameBase<PowerSeriesEndgame<FixedPrecEndgame<tracking::MultiplePrecisionTracker>>, FixedPrecEndgame<tracking::MultiplePrecisionTracker>>;
extern template class EndgameBase<PowerSeriesEndgame<AMPEndgame>, AMPEndgame>;
extern template class EndgameBase<CauchyEndgame<FixedPrecEndgame<tracking::DoublePrecisionTracker>>, FixedPrecEndgame<tracking::DoublePrecisionTracker>>;
extern template class EndgameBase<CauchyEndgame<FixedPrecEndgame<tracking::MultiplePrecisionTracker>>, FixedPrecEndgame<tracking::MultiplePrecisionTracker>>;
extern template class EndgameBase<CauchyEndgame<AMPEndgame>, AMPEndgame>;

extern template class PowerSeriesEndgame<FixedPrecEndgame<tracking::DoublePrecisionTracker>>;
extern template class PowerSeriesEndgame<FixedPrecEndgame<tracking::MultiplePrecisionTracker>>;
extern template class PowerSeriesEndgame<AMPEndgame>;
extern template class CauchyEndgame<FixedPrecEndgame<tracking::DoublePrecisionTracker>>;
extern template class CauchyEndgame<FixedPrecEndgame<tracking::MultiplePrecisionTracker>>;
extern template class CauchyEndgame<AMPEndgame>;

}} // namespaces


