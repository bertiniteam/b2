//This file is part of Bertini 2.
//
//eval_method.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//eval_method.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with eval_method.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file bertini2/system/eval_method.hpp

\brief The EvalMethod enum and its default.

Split out of system.hpp so the evaluation blocks (which carry an eval/deriv method) can see
the full enum definitions without depending on the whole System header — system.hpp includes
the block variant before the System class is defined.
*/

#pragma once

namespace bertini {

	enum class EvalMethod
	{
		FunctionTree, // using virtual methods and recursion
		SLP // using straight line programs
		    // now!  20230714, Eindhoven, Netherlands
	};

	/**
	\brief Gets the default evaluation method for Jacobians.  One might be faster...
	*/
	EvalMethod DefaultEvalMethod();

} // namespace bertini
