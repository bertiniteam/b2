//This file is part of Bertini 2.
//
//node.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//node.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with node.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of notre dame
// Jeb Collins, West Texas A&M


#include "bertini2/function_tree.hpp"

#include <unordered_map>
#include <vector>
#include <mutex>
#include <algorithm>

BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::Variable)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::Differential)

BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::Float)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::Integer)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::Rational)

BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::special_number::Pi)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::special_number::E)

BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::NamedExpression)

BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::SinOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::ArcSinOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::CosOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::ArcCosOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::TanOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::ArcTanOperator)


BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::SumOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::NegateOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::MultOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::PowerOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::IntegerPowerOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::SqrtOperator)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::ExpOperator)




namespace bertini{
namespace node{

	// Default: nothing to simplify -- return this node unchanged (sharing preserved).
	// Operators override to recurse + reassemble through the Simplified* factories.
	std::shared_ptr<Node> Node::Simplified() const
	{
		return std::const_pointer_cast<Node>(shared_from_this());
	}

	// Default: nothing to homogenize (leaves) -- return this node unchanged.  Operators that
	// can carry degree-deficient summands (and their ancestors) override to rebuild functionally.
	std::shared_ptr<Node> Node::Homogenized(VariableGroup const& /*vars*/, std::shared_ptr<Variable> const& /*homvar*/) const
	{
		return std::const_pointer_cast<Node>(shared_from_this());
	}

	// ---- structural hash / equality (predicate layer for hash-consing) ----

	std::size_t Node::Hash() const
	{
		if (!structural_hash_)
			structural_hash_ = HashImpl();
		return *structural_hash_;
	}

	// Default: identity hash (the object address).  Distinct objects hash distinctly; value
	// and operator nodes override HashImpl to be structural.
	std::size_t Node::HashImpl() const
	{
		return std::hash<const void*>{}(this);
	}

	// Default: identity equality.  Value/operator nodes override.
	bool Node::IsSame(Node const& other) const
	{
		return this == &other;
	}


	bool Node::IsPolynomial(std::shared_ptr<Variable> const&v) const
	{
		return Degree(v)>=0;
	}

	bool Node::IsPolynomial(VariableGroup const&v) const
	{
		return Degree(v)>=0;
	}

	Node::Node()
	{ }

	// ---- hash-consing intern table ----
	namespace {
		// Process-global table: structural hash -> live nodes, held weakly so it self-cleans
		// (a node dies when its last external shared_ptr drops; its weak_ptr is pruned on the
		// next touch of that bucket).  Lazy-init function-local statics avoid SIOF.
		std::unordered_map<std::size_t, std::vector<std::weak_ptr<Node>>>& InternBuckets()
		{
			static std::unordered_map<std::size_t, std::vector<std::weak_ptr<Node>>> buckets;
			return buckets;
		}
		std::mutex& InternMutex()
		{
			static std::mutex m;
			return m;
		}
	}

	std::shared_ptr<Node> Intern(std::shared_ptr<Node> const& candidate)
	{
		std::lock_guard<std::mutex> lock(InternMutex());
		auto& bucket = InternBuckets()[candidate->Hash()];

		std::shared_ptr<Node> found;
		// scan for a live, structurally-equal node; prune any expired weak_ptrs as we go
		bucket.erase(
			std::remove_if(bucket.begin(), bucket.end(),
				[&](std::weak_ptr<Node> const& wp) {
					auto sp = wp.lock();
					if (!sp) return true;                          // dead -> prune
					if (!found && sp->IsSame(*candidate)) found = sp;
					return false;
				}),
			bucket.end());

		if (found)
			return found;                                          // hit: discard the candidate
		bucket.push_back(candidate);                               // miss: register and keep
		return candidate;
	}
} // namespace node
} // namespace bertini


