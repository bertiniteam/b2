//This file is part of Bertini 2.0.
//
//symbol.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//symbol.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with symbol.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// symbol.hpp:  Declares the class Symbol.

#ifndef b2Test_Symbol_h
#define b2Test_Symbol_h


#include "function_tree/node.hpp"



namespace  bertini {
	
	/**
	 Node -> Symbol
	
	This class could work as an interface for all non-operators.
	*/
	class Symbol : public virtual Node
	{
		
	public:
		
		virtual ~Symbol() = default;
		
		
	};
	
	
	
	
	
	
	
	
	/**
	Symbols which have names are named symbols.
	*/
	class NamedSymbol : public virtual Symbol
	{
		std::string name_;
	public:
		
		/**
		Get the name of the named symbol
		*/
		std::string name() const
		{
			return name_;
		};
		
		/**
		Get the name of the named symbol
		*/
		void name(const std::string & new_name){name_ = new_name;};
		
		/**
		Default constructor
		*/
		NamedSymbol(){};

		/**
		Parameterized constructor, sets the name of the symbol
		*/
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
	
	
} // re: namespace bertini

#endif
