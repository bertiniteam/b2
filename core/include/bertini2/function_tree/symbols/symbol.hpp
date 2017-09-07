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
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
// Dani Brake
// University of Notre Dame
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
	class Symbol : public virtual Node
	{
		
	public:
		
		virtual ~Symbol() = default;

		unsigned EliminateZeros() override;
		unsigned EliminateOnes() override;
		
	private:
		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Node>(*this);
		}
	};
	
	
	
	
	
	
	
	
	/**
	\brief Symbols which have names are named symbols.
	
	Symbols which have names are named symbols.
	*/
	class NamedSymbol : public virtual Symbol
	{
		std::string name_;
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
	private:

		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned version) {
			ar & boost::serialization::base_object<Symbol>(*this);
			ar & name_;
		}
		
	};
	
} // re: namespace node	
} // re: namespace bertini

#endif
