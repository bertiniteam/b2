//This file is part of Bertini 2.
//
//visitor.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//visitor.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with visitor.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire
//

/**
\file visitor.hpp

\brief Contains the visitor base types.
*/

#ifndef BERTINI_DETAIL_VISITOR_HPP
#define BERTINI_DETAIL_VISITOR_HPP
#include <boost/type_index.hpp>

namespace bertini{

	/**
	\brief A strawman class for implementing the visitor pattern.  

	Deliberately empty, the non-base visitors must derive from it.

	\see Visitor, Observer, AnyObserver, MultiObserver
	*/
	class VisitorBase
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		virtual ~VisitorBase() = default;
	};


	/**
	\brief The first non-trivial visitor class.  

	Derive from this when creating new visitor types.

	\tparam VisitedT The type of object the visitor visits.
	\tparam RetT The type of object to be returned when visiting.  Default is `void`.
	*/
	template<class VisitedT, typename RetT = void>
	class Visitor
	{ BOOST_TYPE_INDEX_REGISTER_CLASS
	public:
		typedef RetT ReturnType;
		virtual ReturnType Visit(VisitedT const &) = 0;
		virtual ~Visitor() = default;
	};

	

}

#endif
