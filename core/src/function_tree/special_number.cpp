//This file is part of Bertini 2.0.
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

//  special_number.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015


#include "function_tree/symbols/special_number.hpp"



namespace bertini {
namespace node{

	std::shared_ptr<Node> Pi()
	{
		return std::make_shared<special_number::Pi>();
	}

	std::shared_ptr<Node> E()
	{
		return std::make_shared<special_number::E>();
	}

	std::shared_ptr<Node> I()
	{
		return std::make_shared<Float>("0.0","1.0");
	}


	std::shared_ptr<Node> Two()
	{
		return std::make_shared<Integer>(2);
	}

	std::shared_ptr<Node> One()
	{
		return std::make_shared<Integer>(1);
	}

	std::shared_ptr<Node> Zero()
	{
		return std::make_shared<Integer>(0);
	}

}
}

