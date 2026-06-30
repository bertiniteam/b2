//This file is part of Bertini 2.
//
//bertini2/blackbox/config.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/blackbox/config.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/blackbox/config.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.
//
// silviana amethyst, university of wisconsin-eau claire

/**
\file bertini2/blackbox/config.hpp 

\brief Configs for the blackbox whatnot
*/

#pragma once

namespace bertini{
	namespace blackbox{


namespace type{
/// \brief Which start system the zero-dim solve should use.
enum class Start{ TotalDegree, RootsOfUnity, MHom, User};
/// \brief Which path tracker the zero-dim solve should use.
enum class Tracker{ FixedDouble, FixedMultiple, Adaptive};
/// \brief Which endgame the zero-dim solve should use.
enum class Endgame{ PowerSeries, Cauchy};
}

// NOTE: the old StorageSelector<StartSystem> (clone-vs-ref per start-system type) is
// gone.  ZeroDim is no longer templated on the start system, so the clone-vs-ref choice
// is made at the construction site by picking the policy instantiation directly: the
// blackbox always clones (CloneGiven + a start-system factory; see switches_zerodim.hpp),
// and user homotopies use RefToGiven through the dedicated UserHomotopy path.

	}
}
