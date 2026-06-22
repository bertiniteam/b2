//This file is part of Bertini 2.
//
//bertini2/function_tree/roots/named_expression.hpp is free software: you can redistribute it
//and/or modify it under the terms of the GNU General Public License as published by the Free
//Software Foundation, either version 3 of the License, or (at your option) any later version.
//
//bertini2/function_tree/roots/named_expression.hpp is distributed in the hope that it will be
//useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License along with
//bertini2/function_tree/roots/named_expression.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire

/**
\file bertini2/function_tree/roots/named_expression.hpp

\brief A user-named expression node: prints as its name, evaluates to its expression.
*/

#ifndef BERTINI_NAMED_EXPRESSION_HPP
#define BERTINI_NAMED_EXPRESSION_HPP

#include "bertini2/function_tree/roots/function.hpp"  // for Handle

namespace bertini {
namespace node{

	/**
	\brief A user-named expression --- the surviving, immutable form of the old Function/Handle.

	`Named(x^2+y^2, "a")` wraps an expression and gives it a name.  It is **immutable**: the
	expression is supplied at construction (there is no SetRoot).  It is hash-consed by
	(expression, name), so `Named(e,"a")` is a distinct node from the bare `e` and from
	`Named(e,"b")`.  It **prints as its name** --- the expansion is revealed elsewhere (the
	System's Describe, which discovers named expressions).  Everything else (eval, differentiate,
	degree, homogenize, precision) forwards to the wrapped expression via Handle.

	Distinct from Variable (a leaf) and from the core's named symbols (Pi, E): those are the other
	two named kinds.  Named expressions are discoverable as their own kind (see Find).
	*/
	class NamedExpression : public Handle
	{
	public:
		BERTINI_DEFAULT_VISITABLE()

		template<typename... Ts>
		static
		std::shared_ptr<NamedExpression> Make(Ts&& ...ts){
			return std::static_pointer_cast<NamedExpression>(Intern(std::shared_ptr<Node>( new NamedExpression(ts...) )));
		}

		/// Prints as just the name (the expansion is shown by the System's Describe).
		void print(std::ostream & target) const override
		{
			target << name();
		}

		// Hash-consed by (expression, name), so distinct names / expressions are distinct nodes.
		std::size_t HashImpl() const override
		{
			std::size_t h = typeid(NamedExpression).hash_code();
			HashCombine(h, EntryNode()->Hash());
			HashCombine(h, std::hash<std::string>{}(name()));
			return h;
		}
		bool IsSame(Node const& other) const override
		{
			auto o = dynamic_cast<NamedExpression const*>(&other);
			// entries are interned, so pointer identity is structural equality
			return o && name() == o->name() && EntryNode().get() == o->EntryNode().get();
		}

		virtual ~NamedExpression() = default;

	protected:
		NamedExpression() = default;

	private:
		NamedExpression(std::shared_ptr<Node> const& entry, std::string const& name)
			: Handle(entry, name)
		{}

		friend class boost::serialization::access;

		template <typename Archive>
		void serialize(Archive& ar, const unsigned /*version*/) {
			ar & boost::serialization::base_object<Handle>(*this);
		}
	};

	/// Give an expression a name: `Named(x*x + y*y, "a")`.  The result prints as `a` and evaluates
	/// to the expression; it is hash-consed by (expression, name).
	inline std::shared_ptr<NamedExpression> Named(std::shared_ptr<Node> const& expr, std::string const& name)
	{
		return NamedExpression::Make(expr, name);
	}

} // re: namespace node
} // re: namespace bertini

#endif
