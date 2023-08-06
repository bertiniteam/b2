//This file is part of Bertini 2.
//
//visitable.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//visitable.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with visitable.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
\file visitable.hpp

\brief Contains the visitable base types
*/

#ifndef BERTINI_DETAIL_VISITABLE_HPP
#define BERTINI_DETAIL_VISITABLE_HPP

#include "bertini2/detail/visitor.hpp"
#include "bertini2/detail/events.hpp"

namespace bertini{

	namespace policy{

		/**
		\brief A policy for what to do when a visited type is encoutered by a visitor it doesn't know.

		Defines the default behaviour when an unknown visitor is encountered during visitation.  

		Assumes the return type is default constructible.
		*/
		template<class VisitedT, typename RetT>
		struct DefaultConstruct
		{
			static RetT OnUnknownVisitor(VisitedT&a, VisitorBase&) // VisitorBase is inherited from using CRTP -- it's not a template parameter to this function, so its typeid name is useless.
			{
				std::cout << "unknown visitor: " << a << " of type " << typeid(VisitedT).name() << ".  Make sure you've added it in all three places!  two in the class definition (type inheritance listing as visiting, and virtual function declaration), and one in cpp source (function definition)." << std::endl;
				return VisitedT();
			}
		};

		template<class VisitedT, typename RetT>
		struct RaiseExceptionWithTypeNamesInMessage
		{
			static RetT OnUnknownVisitor(VisitedT&a, VisitorBase&)
			{
				std::stringstream err_msg;
				err_msg << "unknown visitor: " << a << " of type " << typeid(VisitedT).name() << ".  Make sure you've added it in all three places!  two in the class definition (type inheritance listing as visiting, and virtual function declaration), and one in cpp source (function definition)." << std::endl;
				throw std::runtime_error(err_msg.str());
			}
		};
		/**
		The default policy for what to do when a visitable is visited by an unknown visitor.
		*/
		template<class VisitedT, typename RetT>
		using DefaultCatchAll = RaiseExceptionWithTypeNamesInMessage<VisitedT, RetT>;
	} // namespace policy

	


	/**
	\brief The base class for visitable types.  

	Implemented based on Alexandrescu, 2001, and Hythem Sidky's SAPHRON package, with his permission.

	\see BERTINI_DEFAULT_VISITABLE
	*/
	template< typename RetT = void, template<class,class> class CatchAll = policy::DefaultCatchAll>
	class VisitableBase
	{
	public:
		typedef RetT VisitReturnType;
		virtual ~VisitableBase() = default;
		virtual VisitReturnType Accept(VisitorBase&) = 0; // the implementation will either be provided by a macro, or by hand, for each class which is visitable.

	protected:

		/**
		\brief Abstract method for how to accept a visitor.  

		This function simply Forwards to the Visit method of the visitor, if the visitor's type is known.  If known, invokes the behaviour defined by the CatchAll template parameter.

		\tparam T The type of the visited object.  Should be inferred by the compiler.

		\param visited The visited visitable object.
		\param guest The visiting object.
		*/
		template<typename T>
		static
		VisitReturnType AcceptBase(T& visited, VisitorBase& guest)
		{
			if (auto p = dynamic_cast<Visitor<T>*>(&guest))
				return p->Visit(visited);
			else
				return CatchAll<T,RetT>::OnUnknownVisitor(visited, guest);
		}
	};

	
	/** 
	\brief macro for classes which want default Accept implementation, having nothing fancy to do when accepting.
	*/
	#define BERTINI_DEFAULT_VISITABLE() \
		virtual VisitReturnType Accept(VisitorBase& guest) override \
		{ return AcceptBase(*this, guest); }

} // namespace bertini


#endif

