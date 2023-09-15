//This file is part of Bertini 2.
//
//special_number.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//special_number.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with special_number.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015-2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire

#include "bertini2/function_tree/symbols/special_number.hpp"



namespace bertini {
namespace node{

	std::shared_ptr<Node> Pi()
	{
		return special_number::Pi::Make();
	}

	std::shared_ptr<Node> E()
	{
		return special_number::E::Make();
	}

	std::shared_ptr<Node> I()
	{
		return MakeFloat(0,1);
	}


	std::shared_ptr<Node> Two()
	{
		return MakeInteger(2);
	}

	std::shared_ptr<Node> One()
	{
		return MakeInteger(1);
	}

	std::shared_ptr<Node> Zero()
	{
		return MakeInteger(0);
	}

}
}

