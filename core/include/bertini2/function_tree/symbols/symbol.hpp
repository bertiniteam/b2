//This file is part of Bertini 2.
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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
// silviana amethyst, university of wisconsin-eau claire
//
//  Created by Collins, James B. on 4/30/15.
//
//
// symbol.hpp:  Declares the class Symbol.

/**
\file symbol.hpp

\brief Defines the abstract Symbol and NamedSymbol classes,

*/

#ifndef BERTINI_FUNCTION_TREE_SYMBOL_HPP
#define BERTINI_FUNCTION_TREE_SYMBOL_HPP


#include "bertini2/function_tree/node.hpp"



namespace  bertini {
namespace node {
	/**
	 \brief Abstract symbol class.

	This class is an interface for all non-operators.
	*/
	class Symbol : public Node
	{
		
	public:
		
		virtual ~Symbol() = default;

	private:
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<Node>(*this);
		}
	};
	
	
	
	
	
	
	
	
	/**
	\brief Capability class for things which have a name.

	Deliberately NOT a Node: classes that need a name alongside a different
	primary base (e.g. special_number::Pi, which is a Number) inherit this
	without creating a diamond in the Node hierarchy.
	*/
	class Named
	{
	public:
		const std::string& name() const
		{
			return name_;
		}

		void name(const std::string& new_name)
		{
			name_ = new_name;
		}

	protected:
		~Named() = default;  // not polymorphic; never delete through Named*

		Named() = default;

		explicit Named(std::string new_name) : name_(std::move(new_name))
		{}

		std::string name_;

	private:
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & name_;
		}
	};




	/**
	\brief Symbols which have names are named symbols.

	Symbols which have names are named symbols.
	*/
	class NamedSymbol : public Symbol
	{

	public:
		
		/**
		Get the name of the named symbol
		*/
		const std::string & name() const;
		
		/**
		Get the name of the named symbol
		*/
		void name(const std::string & new_name);
		
		
		/**
		Parameterized constructor, sets the name of the symbol
		*/
		NamedSymbol(const std::string & new_name);
		
		
		void print(std::ostream& target) const override;
		
		virtual ~NamedSymbol() = default;
		
	protected:
		NamedSymbol() = default;

		std::string name_;

	private:

		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<Symbol>(*this);
			ar & name_;
		}
		
	};
	
} // re: namespace node	
} // re: namespace bertini

#endif
