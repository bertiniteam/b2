//This file is part of Bertini 2.
//
//canonical.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//canonical.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with canonical.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

#include "bertini2/function_tree.hpp"
#include "bertini2/function_tree/canonical.hpp"
#include "bertini2/function_tree/gather.hpp"

#include <set>
#include <string>
#include <sstream>
#include <numeric>
#include <algorithm>

namespace bertini {
namespace node {

namespace {
	// session-global canonicalization settings.  Canonicalization is ON by default: every Sum/Mult normalizes operand
	// order, so structurally-equal expressions dedup to one interned node.
	MonomialOrder& TheOrder()        { static MonomialOrder o = MonomialOrder::GrevLex; return o; }
	bool&          TheEnabledFlag()  { static bool on = true; return on; }
	bool&          ThePowerFoldFlag(){ static bool on = true; return on; }

	bool AllNonNegative(std::vector<int> const& v)
	{
		for (int e : v) if (e < 0) return false;
		return true;
	}

	// a constant factor: a polynomial monomial of total degree zero (no variables).
	bool IsConstantKey(std::vector<int> const& v)
	{
		return AllNonNegative(v) && std::accumulate(v.begin(), v.end(), 0) == 0;
	}

	// Does monomial 'a' sort BEFORE monomial 'b' (i.e. is it the "greater" / leading one,
	// since we sort descending)?  Polynomial monomials precede non-polynomial operands.
	bool MonomialGreater(std::vector<int> const& a, std::vector<int> const& b, MonomialOrder order)
	{
		const bool aPoly = AllNonNegative(a);
		const bool bPoly = AllNonNegative(b);
		if (aPoly != bPoly) return aPoly;     // polynomial terms come first
		if (!aPoly)         return false;     // both non-polynomial: a tie (printed-form breaks it)

		const std::size_t n = a.size();
		switch (order)
		{
			case MonomialOrder::Lex:
				for (std::size_t i = 0; i < n; ++i)
					if (a[i] != b[i]) return a[i] > b[i];
				return false;

			case MonomialOrder::RevLex:
				for (std::size_t i = n; i-- > 0; )
					if (a[i] != b[i]) return a[i] < b[i];
				return false;

			case MonomialOrder::GrevLex:
			{
				const int da = std::accumulate(a.begin(), a.end(), 0);
				const int db = std::accumulate(b.begin(), b.end(), 0);
				if (da != db) return da > db;             // graded: higher total degree first
				for (std::size_t i = n; i-- > 0; )         // then reverse-lex on ties
					if (a[i] != b[i]) return a[i] < b[i];
				return false;
			}
		}
		return false;
	}

	std::string PrintOf(std::shared_ptr<Node> const& n)
	{
		std::ostringstream oss;
		n->print(oss);
		return oss.str();
	}

	// A product operand decomposed into (base, integer exponent): x -> (x,1), x^k -> (x,k).
	std::pair<std::shared_ptr<Node>, int> BaseAndExponent(std::shared_ptr<Node> const& n)
	{
		if (auto ip = std::dynamic_pointer_cast<IntegerPowerOperator>(n))
			return { ip->Operand(), ip->exponent() };
		return { n, 1 };
	}

	// A base whose repeated occurrences we collapse into a single power.  Numeric constants are
	// left alone (2*2 stays 2*2 -- constant folding is a separate pass, not 2^2); everything else
	// (variables and compound subexpressions like (x+y) or sin(x)) folds -- that is the CSE win.
	bool IsFoldableBase(std::shared_ptr<Node> const& n)
	{
		return !std::dynamic_pointer_cast<Number>(n);
	}

	// Flatten a product's operand list: splice any operand that is itself a MultOperator into this
	// list, composing the mult/div flag (a child factor's flag cf under a parent slot flag pf
	// becomes pf==cf, i.e. div-of-div = mult).  Fully recursive; single-operand MultOperator
	// wrappers are unwrapped here too.
	void FlattenProduct(std::vector<std::shared_ptr<Node>>& operands, std::vector<bool>& flags)
	{
		std::vector<std::shared_ptr<Node>> out_ops;
		std::vector<bool>                  out_flags;
		std::vector<std::pair<std::shared_ptr<Node>, bool>> stack;   // pre-order traversal
		stack.reserve(operands.size());
		for (std::size_t i = operands.size(); i-- > 0; )
			stack.emplace_back(operands[i], flags[i]);

		while (!stack.empty())
		{
			auto node   = stack.back().first;
			const bool f = stack.back().second;
			stack.pop_back();

			// Splice a nested product only when it is a MULTIPLICAND.  A divisor sub-product like
			// x/(y*z) must stay atomic: splicing it would distribute the division into x/y/z (two
			// divides) instead of one multiply plus one divide -- a pessimization, division being the
			// costliest op.  Multiplicand nesting (the common y*x*x case) still fully flattens.
			auto mo = std::dynamic_pointer_cast<MultOperator>(node);
			if (mo && f)
			{
				auto const& kids   = mo->Operands();
				auto const& kflags = mo->GetMultOrDiv();
				for (std::size_t j = kids.size(); j-- > 0; )
					stack.emplace_back(kids[j], kflags[j]);   // f is true here, so f==kflags[j] == kflags[j]
			}
			else
			{
				out_ops.push_back(node);
				out_flags.push_back(f);
			}
		}
		operands.swap(out_ops);
		flags.swap(out_flags);
	}

	// Fold like factors of a (flattened) product: group operands by base, sum the signed exponents
	// (multiplicand +, divisor -), and emit one factor per base -- x*x -> x^2, x*x/x -> x, x/x -> 1.
	// Numeric-constant operands are passed through untouched, each its own group.
	void FoldLikeFactors(std::vector<std::shared_ptr<Node>>& operands, std::vector<bool>& flags)
	{
		struct Group {
			std::shared_ptr<Node> base;      // foldable base (null for a passthrough constant)
			std::string           key;       // printed form of base, for matching
			int                   exp;       // net signed exponent
			bool                  foldable;
			std::shared_ptr<Node> passthru;  // the original operand, for constants
			bool                  passflag;
		};

		std::vector<Group> groups;
		for (std::size_t i = 0; i < operands.size(); ++i)
		{
			auto be = BaseAndExponent(operands[i]);
			if (IsFoldableBase(be.first))
			{
				const int signed_exp = flags[i] ? be.second : -be.second;
				std::string k = PrintOf(be.first);
				bool merged = false;
				for (auto& g : groups)
					if (g.foldable && g.key == k) { g.exp += signed_exp; merged = true; break; }
				if (!merged)
					groups.push_back(Group{ be.first, std::move(k), signed_exp, true, nullptr, false });
			}
			else
				groups.push_back(Group{ nullptr, std::string{}, 0, false, operands[i], flags[i] });
		}

		std::vector<std::shared_ptr<Node>> out_ops;
		std::vector<bool>                  out_flags;
		for (auto& g : groups)
		{
			if (!g.foldable)
			{
				out_ops.push_back(g.passthru);
				out_flags.push_back(g.passflag);
				continue;
			}
			if (g.exp == 0)                                 // net factor of 1 -- drops out
				continue;
			const int a = g.exp < 0 ? -g.exp : g.exp;
			std::shared_ptr<Node> node = (a == 1)
				? g.base
				: std::static_pointer_cast<Node>(IntegerPowerOperator::Make(g.base, a));
			out_ops.push_back(node);
			out_flags.push_back(g.exp > 0);
		}
		if (out_ops.empty())                                // everything cancelled -> the product is 1
		{
			out_ops.push_back(std::static_pointer_cast<Node>(Integer::Make(1)));
			out_flags.push_back(true);
		}
		operands.swap(out_ops);
		flags.swap(out_flags);
	}
} // anon namespace

MonomialOrder CurrentMonomialOrder()           { return TheOrder(); }
void          SetMonomialOrder(MonomialOrder o){ TheOrder() = o; }
bool          CanonicalizeByDefault()          { return TheEnabledFlag(); }
void          SetCanonicalizeByDefault(bool on){ TheEnabledFlag() = on; }
bool          PowerFoldByDefault()             { return ThePowerFoldFlag(); }
void          SetPowerFoldByDefault(bool on)   { ThePowerFoldFlag() = on; }

void CanonicalizeNaryOperands(std::vector<std::shared_ptr<Node>>& operands,
                              std::vector<bool>& flags,
                              bool multiplicands_first)
{
	if (!CanonicalizeByDefault() || operands.size() < 2)
		return;

	// A product is normalized to a flat, power-folded form before sorting: flatten nested
	// MultOperators into one factor list, then collapse repeated bases into IntegerPowers
	// (x*x -> x^2, y*x*x -> x^2*y, x*x/x -> x).  This is what makes a squared factor the SAME
	// interned IntegerPower node the differentiator emits, so the SLP computes it once and shares
	// it between the function and its Jacobian.  Sums keep their existing sort-only behavior.
	if (multiplicands_first && PowerFoldByDefault())
	{
		FlattenProduct(operands, flags);
		FoldLikeFactors(operands, flags);
		if (operands.size() < 2)   // fold may have collapsed the product to a single factor
			return;
	}

	// global variable order: the union of the operands' variables, alphabetical by name
	// (GatherVariables already sorts by name; variables are canonical-by-name post-3b).
	VariableGroup vars;
	std::set<std::string> seen;
	for (auto const& op : operands)
		for (auto const& v : GatherVariables(op))
			if (seen.insert(v->name()).second)
				vars.push_back(v);
	std::sort(vars.begin(), vars.end(),
	          [](std::shared_ptr<Variable> const& a, std::shared_ptr<Variable> const& b)
	          { return a->name() < b->name(); });

	const std::size_t N = operands.size();
	std::vector<std::vector<int>> keys(N);
	std::vector<std::string>      prints(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		keys[i]   = operands[i]->MultiDegree(vars);
		prints[i] = PrintOf(operands[i]);
	}

	const MonomialOrder order = CurrentMonomialOrder();
	std::vector<std::size_t> perm(N);
	std::iota(perm.begin(), perm.end(), std::size_t{0});
	std::stable_sort(perm.begin(), perm.end(),
		[&](std::size_t i, std::size_t j) -> bool
		{
			// for a Mult, multiplicands (flag true) sort ahead of divisors (flag false), so a
			// divisor is never first; the +/- sign of a Sum term does not affect ordering.
			if (multiplicands_first && flags[i] != flags[j])
				return flags[i];
			// within a product, a constant coefficient sorts first, so monomials read the
			// conventional way: "3*x^2", not "x^2*3".  (Sums keep degree order: the constant
			// term stays last, e.g. "x^2+2*x*y-1".)
			if (multiplicands_first)
			{
				const bool ci = IsConstantKey(keys[i]);
				const bool cj = IsConstantKey(keys[j]);
				if (ci != cj) return ci;
			}
			if (MonomialGreater(keys[i], keys[j], order)) return true;
			if (MonomialGreater(keys[j], keys[i], order)) return false;
			return prints[i] < prints[j];   // deterministic, content-based tie-break
		});

	std::vector<std::shared_ptr<Node>> new_operands(N);
	std::vector<bool>                  new_flags(N);
	for (std::size_t i = 0; i < N; ++i)
	{
		new_operands[i] = operands[perm[i]];
		new_flags[i]    = flags[perm[i]];
	}
	operands.swap(new_operands);
	flags.swap(new_flags);
}

} // namespace node
} // namespace bertini
