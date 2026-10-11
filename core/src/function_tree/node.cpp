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

#include "bertini2/function_tree.hpp"
#include "bertini2/function_tree/canonical_encoding.hpp"
#include "bertini2/detail/sha256.hpp"

#include <memory>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <algorithm>

BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::Variable)
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::Differential)

BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::Complex)
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
BOOST_CLASS_EXPORT_IMPLEMENT(bertini::node::LogOperator)




namespace bertini{
namespace node{

    // Default: nothing to simplify -- return this node unchanged (sharing preserved).
    // Operators override to recurse + reassemble through the Simplified* factories.
    std::shared_ptr<Node> Node::Simplified() const
    {
        return std::const_pointer_cast<Node>(shared_from_this());
    }

    // Default: nothing to substitute (leaves other than Variable, and Differentials) -- return
    // this node unchanged.  Variable overrides to match; operators/NamedExpression override to
    // recurse + reassemble through the Simplified* factories.
    std::shared_ptr<Node> Node::Subs(SubstitutionMap const& /*substitutions*/) const
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

    // ---- differentiation, each distinct node once (#505) ----

    // The derivatives kept by the outermost open DerivativeScope on a thread, keyed by
    // (node, variable).  Each entry holds its node and variable as well as the derivative: a key
    // is an address, and a node that died while the scope was open could have its address reused
    // by another node, which would then be answered with the dead node's derivative.
    struct DerivativeScope::Derivatives
    {
        struct Key
        {
            Node const* node;
            Variable const* variable;
            bool operator==(Key const& other) const { return node == other.node && variable == other.variable; }
        };
        struct KeyHash
        {
            std::size_t operator()(Key const& k) const
            {
                std::size_t seed = std::hash<Node const*>{}(k.node);
                seed ^= std::hash<Variable const*>{}(k.variable) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
                return seed;
            }
        };
        struct Entry
        {
            std::shared_ptr<const Node> node;
            std::shared_ptr<Variable> variable;
            std::shared_ptr<Node> derivative;
        };
        std::unordered_map<Key, Entry, KeyHash> entries;
    };

    namespace {
        thread_local DerivativeScope::Derivatives* open_derivatives = nullptr;
    }

    DerivativeScope::DerivativeScope()
    {
        if (!open_derivatives)
        {
            owned_ = std::make_unique<Derivatives>();
            open_derivatives = owned_.get();
        }
    }

    DerivativeScope::~DerivativeScope()
    {
        if (owned_)
            open_derivatives = nullptr;
    }

    std::shared_ptr<Node> Node::Differentiate(std::shared_ptr<Variable> const& v) const
    {
        DerivativeScope scope;
        auto& entries = open_derivatives->entries;
        DerivativeScope::Derivatives::Key const key{this, v.get()};
        auto found = entries.find(key);
        if (found != entries.end())
            return found->second.derivative;
        auto derivative = DifferentiateImpl(v);
        // the rule differentiated other nodes meanwhile and may have rehashed the table, so
        // insert afresh rather than through anything found before it ran
        entries.emplace(key, DerivativeScope::Derivatives::Entry{shared_from_this(), v, derivative});
        return derivative;
    }

    // Convenience: repeated single-variable differentiation.  Folds the one-variable
    // Differentiate `count` times; count==0 is the identity (the node itself).
    std::shared_ptr<Node> Node::Differentiate(std::shared_ptr<Variable> const& v, unsigned count) const
    {
        if (count==0)
            return std::const_pointer_cast<Node>(shared_from_this());
        auto result = Differentiate(v);
        for (unsigned i = 1; i < count; ++i)
            result = result->Differentiate(v);
        return result;
    }

    // Convenience: sequential differentiation wrt each variable in the group (mixed partials).
    // Folds the one-variable Differentiate over `vars` in order; empty is the identity.
    std::shared_ptr<Node> Node::Differentiate(VariableGroup const& vars) const
    {
        if (vars.empty())
            return std::const_pointer_cast<Node>(shared_from_this());
        auto result = Differentiate(vars.front());
        for (size_t i = 1; i < vars.size(); ++i)
            result = result->Differentiate(vars[i]);
        return result;
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

    // ---- the degree record: variables, degrees and multidegrees, once per node (#505) ----

    struct Node::DegreeRecord
    {
        VariableGroup variables;               ///< distinct, sorted by name (GatherVariables)
        VariableGroup mentioned;               ///< variables and differentials' variables, sorted by address
        std::vector<int> degree;               ///< Degree(v) for each mentioned v
        std::vector<int> multidegree;          ///< the MultiDegree entry for each mentioned v
        int total_degree = 0;                  ///< Degree(nullptr)
        int unmentioned_degree = 0;            ///< Degree(v) for any v not mentioned
        int unmentioned_multidegree = 0;       ///< the MultiDegree entry for any v not mentioned

        // the index of v among the mentioned variables, or -1
        long IndexOf(Variable const* v) const
        {
            auto it = std::lower_bound(mentioned.begin(), mentioned.end(), v,
                                       [](std::shared_ptr<Variable> const& a, Variable const* b) { return a.get() < b; });
            return (it != mentioned.end() && it->get() == v) ? it - mentioned.begin() : -1;
        }
    };

    // A variable no expression contains: asking a node about it asks about every variable the
    // node does not mention.  Made with the default constructor and never interned, so no
    // variable a user makes, whatever its name, is this one.
    std::shared_ptr<Variable> const& Node::Unmentioned()
    {
        static std::shared_ptr<Variable> const unmentioned(new Variable());
        return unmentioned;
    }

    namespace {

        // Each operand of n, by the public child accessors (the same descent as Find).
        template<typename F>
        void ForEachOperand(Node const& n, F&& f)
        {
            if (auto h = dynamic_cast<NamedExpression const*>(&n))
                f(h->EntryNode());
            else if (auto nary = dynamic_cast<NaryOperator const*>(&n))
                for (auto const& child : nary->Operands())
                    f(child);
            else if (auto u = dynamic_cast<UnaryOperator const*>(&n))
                f(u->Operand());
            else if (auto p = dynamic_cast<PowerOperator const*>(&n))
            {
                f(p->GetBase());
                f(p->GetExponent());
            }
        }

        void SortAndDedupe(VariableGroup& vars, bool by_name)
        {
            if (by_name)
                std::sort(vars.begin(), vars.end(),
                          [](std::shared_ptr<Variable> const& a, std::shared_ptr<Variable> const& b) { return a->name() < b->name(); });
            else
                std::sort(vars.begin(), vars.end(),
                          [](std::shared_ptr<Variable> const& a, std::shared_ptr<Variable> const& b) { return a.get() < b.get(); });
            vars.erase(std::unique(vars.begin(), vars.end(),
                                   [](std::shared_ptr<Variable> const& a, std::shared_ptr<Variable> const& b) { return a.get() == b.get(); }),
                       vars.end());
        }
    }

    Node::DegreeRecord const& Node::Record() const
    {
        if (auto made = std::atomic_load(&degree_record_))
            return *made;

        auto record = std::make_shared<DegreeRecord>();
        ForEachOperand(*this, [&](auto const& operand)
        {
            if (!operand)
                return;
            auto const& r = operand->Record();
            record->variables.insert(record->variables.end(), r.variables.begin(), r.variables.end());
            record->mentioned.insert(record->mentioned.end(), r.mentioned.begin(), r.mentioned.end());
        });
        if (dynamic_cast<Variable const*>(this))
        {
            auto as_variable = std::const_pointer_cast<Variable>(std::static_pointer_cast<Variable const>(shared_from_this()));
            record->variables.push_back(as_variable);
            record->mentioned.push_back(as_variable);
        }
        if (auto d = dynamic_cast<Differential const*>(this))
            if (d->GetVariable())
                record->mentioned.push_back(std::const_pointer_cast<Variable>(d->GetVariable()));
        SortAndDedupe(record->variables, true);
        SortAndDedupe(record->mentioned, false);

        // The rules ask the operands' Degree and MultiDegree, which their records answer, and
        // never this node's own (DegreeImpl's contract), so nothing here recurses into this.
        record->degree.reserve(record->mentioned.size());
        for (auto const& v : record->mentioned)
            record->degree.push_back(this->DegreeImpl(v));
        record->total_degree = this->DegreeImpl(nullptr);
        record->unmentioned_degree = this->DegreeImpl(Unmentioned());
        record->multidegree = this->MultiDegreeImpl(record->mentioned);
        record->unmentioned_multidegree = this->MultiDegreeImpl(VariableGroup{Unmentioned()}).front();

        // publish once: a thread that lost the race uses the winner's, so a reference into the
        // record stays valid for the node's life
        std::shared_ptr<const DegreeRecord> expected;
        std::shared_ptr<const DegreeRecord> mine = record;
        if (!std::atomic_compare_exchange_strong(&degree_record_, &expected, mine))
            return *expected;
        return *mine;
    }

    VariableGroup const& Node::Variables() const
    {
        return Record().variables;
    }

    int Node::Degree(std::shared_ptr<Variable> const& v) const
    {
        auto const& record = Record();
        if (!v)
            return record.total_degree;
        auto index = record.IndexOf(v.get());
        return index < 0 ? record.unmentioned_degree : record.degree[index];
    }

    std::vector<int> Node::MultiDegree(VariableGroup const& vars) const
    {
        auto const& record = Record();
        std::vector<int> degrees;
        degrees.reserve(vars.size());
        for (auto const& v : vars)
        {
            auto index = record.IndexOf(v.get());
            degrees.push_back(index < 0 ? record.unmentioned_multidegree : record.multidegree[index]);
        }
        return degrees;
    }

    std::string const& Node::ContentKey() const
    {
        if (auto made = std::atomic_load(&content_key_))
            return *made;
        // the operands' keys are made first, by the shallow encoding asking for them, so a key
        // costs O(operands): a Merkle hash of the canonical encoding
        std::shared_ptr<const std::string> mine = std::make_shared<const std::string>(
            detail::Sha256(ShallowCanonicalEncoding(std::const_pointer_cast<Node>(shared_from_this()))).Hex());
        std::shared_ptr<const std::string> expected;
        if (!std::atomic_compare_exchange_strong(&content_key_, &expected, mine))
            return *expected;
        return *mine;
    }

} // namespace node
} // namespace bertini
