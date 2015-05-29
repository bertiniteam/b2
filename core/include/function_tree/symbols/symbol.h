//This file is part of Bertini 2.0.
//
//Foobar is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//Foobar is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// symbol.h:  Declares the class Symbol.

#ifndef b2Test_Symbol_h
#define b2Test_Symbol_h


#include "function_tree/node.h"





// Node -> Symbol
// Description: This class could work as an interface for all non-operators.
// NOTE: Currently this does nothing!!!!
//
// TODO:(JBC) This defined the FreshEval for all symbols.
class Symbol : public virtual Node
{

public:

	virtual ~Symbol() = default;
	
	
};









class NamedSymbol : public virtual Symbol
{
	std::string name_;
public:
	
	std::string name() const
	{
		return name_;
	};
	
	void name(const std::string & new_name){name_ = new_name;};
	
	NamedSymbol(){};
	NamedSymbol(const std::string & new_name) : name_(new_name){};
	
	
	virtual void print(std::ostream& target) const override
	{
		target << name();
	}
	
	virtual ~NamedSymbol() = default;
	
	
	virtual void Reset() override
	{
		Node::ResetStoredValues();
	}
	
	
};

#endif
