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
	// session-global canonicalization settings (default OFF -- enabling it by default is a
	// deliberate, churn-bearing follow-up step).
	MonomialOrder& TheOrder()        { static MonomialOrder o = MonomialOrder::GrevLex; return o; }
	bool&          TheEnabledFlag()  { static bool on = false; return on; }

	bool AllNonNegative(std::vector<int> const& v)
	{
		for (int e : v) if (e < 0) return false;
		return true;
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
} // anon namespace

MonomialOrder CurrentMonomialOrder()           { return TheOrder(); }
void          SetMonomialOrder(MonomialOrder o){ TheOrder() = o; }
bool          CanonicalizeByDefault()          { return TheEnabledFlag(); }
void          SetCanonicalizeByDefault(bool on){ TheEnabledFlag() = on; }

void CanonicalizeNaryOperands(std::vector<std::shared_ptr<Node>>& operands,
                              std::vector<bool>& flags,
                              bool multiplicands_first)
{
	if (!CanonicalizeByDefault() || operands.size() < 2)
		return;

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
