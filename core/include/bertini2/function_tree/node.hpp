//This file is part of Bertini 2.
//
//node.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//node.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with node.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

/**
\file node.hpp

\brief Defines the abstract Node base class.

*/

#ifndef BERTINI_NODE_BASE_HPP
#define BERTINI_NODE_BASE_HPP




#include <iostream>
#include <string>
#include <map>
#include <tuple>
#include <optional>
#include <typeinfo>

#include <boost/type_index.hpp>

#include "bertini2/num_traits.hpp"
#include "bertini2/detail/visitable.hpp"
#include "bertini2/function_tree/forward_declares.hpp"
#include "bertini2/function_tree/print_dialect.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/nvp.hpp>





namespace bertini {

    namespace node{
        class Variable;
    }

/// An ordered group of variables (the unit of homogenization / variable grouping).
using VariableGroup = std::vector< std::shared_ptr<node::Variable> >;

/// A substitution request for Node::Subs: variable name -> the node to put in its place.
using SubstitutionMap = std::map< std::string, std::shared_ptr<node::Node> >;



/// \brief How a variable group participates in homogenization.
enum class VariableGroupType
{
    Homogeneous,  ///< A projective (homogeneous) variable group.
    Affine,       ///< An affine variable group (to be homogenized).
    Ungrouped     ///< Variables not assigned to any group.
};


namespace node{

/**
\brief Operator precedence classes, used to decide parenthesization when printing.

Higher binds tighter.  A child is wrapped in parentheses only when its
precedence is too low for the position it is printed in; leaves and
self-delimiting nodes (function calls like sin(...), complex pairs) are
atoms and never wrapped.  Negative literal constants report PrecNegate so
they parenthesize exactly where a Negate node would.
*/
enum PrintPrecedence : unsigned {
    PrecSum = 10,
    PrecNegate = 15,
    PrecMult = 20,
    PrecPower = 30,
    PrecAtom = 100
};

/**
An interface for all nodes in a function tree, and for a function object as well.  Almost all
 methods that will be called on a node must be declared in this class.  The main evaluation method is
 defined in this class.

 Bertini function trees are created using std::shared_ptr's to Node base class, generally.

 \brief Abstract base class for the Bertini hybrid-precision (double-multiple) expression tree.
 */
class Node : public VisitableBase<>, public std::enable_shared_from_this<Node>
{
public:

    virtual ~Node() = default;

    ///////// PUBLIC PURE METHODS /////////////////



    /**
    \brief Functionally simplify this tree, returning a NEW simplified tree.

    Non-mutating successor to the in-place EliminateZeros / EliminateOnes / ReduceDepth
    machinery (ADR-0011, issue 251).  The input tree is never modified; a fresh simplified
    tree is returned, built through the SimplifiedSum / SimplifiedMult / SimplifiedNegate
    factories so that literal zeros/ones vanish and exact constants fold.  The default
    (leaves, and any node with nothing to simplify) returns the node unchanged -- structural
    sharing preserved.

    \return A simplified tree (possibly the same node, if nothing simplified).
    */
    virtual std::shared_ptr<Node> Simplified() const;

    /**
    \brief Symbolically substitute variables, returning a NEW tree.

    Replaces each Variable leaf named in `substitutions` with its associated node, functionally
    (the input tree is never modified).  Substitution is simultaneous and single-pass: a leaf is
    replaced at most once and the replacement is not itself descended into, so `{x:y, y:z}` maps
    `x+y` to `y+z` (no cascade) and the order of the map is irrelevant.  A variable absent from
    the map is left unchanged.  The result is rebuilt through the simplifying factories (like
    differentiation), so constants fold.  The default (leaves with nothing to replace) returns the
    node unchanged -- structural sharing preserved.

    Only Variable leaves are matched (variable -> node substitution); substituting into a Jacobian
    tree containing Differentials is not supported.

    \param substitutions A map from variable name to the node to put in that variable's place.
    \return The substituted tree (possibly the same node, if nothing matched).
    */
    virtual std::shared_ptr<Node> Subs(SubstitutionMap const& substitutions) const;

    /**
    Virtual method for printing Nodes to arbitrary output streams.
    */
    virtual void print(std::ostream& target) const = 0;

    /**
    \brief The printing precedence of this node, deciding parenthesization.

    Defaults to PrecAtom: leaves and self-delimiting nodes are never wrapped.
    Operator nodes override this; printing parents wrap a child only when its
    precedence is too low for the position it occupies.
    */
    virtual unsigned Precedence() const
    {
        return PrecAtom;
    }

    /**
    \brief Is this node the literal constant 0?

    Exact check on stored literal values (Integer/Complex/Rational override);
    false for everything else, including non-literal expressions that happen
    to evaluate to zero.  Used by differentiation to build already-simplified
    trees without evaluating anything.
    */
    virtual bool IsLiteralZero() const
    {
        return false;
    }

    /**
    \brief Is this node the literal constant 1?

    \see IsLiteralZero
    */
    virtual bool IsLiteralOne() const
    {
        return false;
    }

    /**
    \brief Memoized structural hash of this node.

    Order-sensitive: equal structures hash equal, where children are folded in by their own
    Hash() and operand order / signs / exponents participate.  The default (leaves and any
    node not overriding HashImpl) is node identity (the object's address), so distinct objects
    hash distinctly; value nodes (Integer/Rational/Complex) and operators override to be
    structural.  Memoized -- valid because nodes are immutable post-construction -- and
    deliberately independent of the mutable working precision (stable across precision()).

    This is the predicate layer for the hash-consing intern table; nothing wires it into
    construction yet.
    */
    std::size_t Hash() const;

    /**
    \brief Opaque, lazily-populated cache of a compiled evaluator for this expression.

    The system layer (EvalExpression / f.eval) compiles a StraightLineProgram for ([this],
    canonical-by-name variable order) on first use and stashes it here, so repeated evaluations
    of the same (hash-consed, immutable) expression reuse it --- like the memoized Hash().  The
    type is erased because the SLP lives in the system layer, above function_tree; the compiled
    program holds no node pointers (its constants are value recipes, ADR-0027), so
    there is no reference cycle.  Transient: not serialized; a clone recompiles on demand.
    */
    std::shared_ptr<const void> EvalProgram() const { return eval_program_; }
    /// \brief Stash a compiled evaluator for this expression (type-erased; see EvalProgram()).
    void SetEvalProgram(std::shared_ptr<const void> p) const { eval_program_ = p; }

    /**
    \brief Order-sensitive structural equality, shallow given interned children.

    Two nodes are the same iff they have the same dynamic type, the same operator payload
    (signs / mult-or-div flags / integer exponent / literal value), and the same operands
    **by pointer** (operands are not recursed -- in the hash-consed world children are already
    canonical, so pointer-equality is structural equality).  The default is node identity;
    value/operator nodes override.  Consistent with Hash(): IsSame(a,b) implies a.Hash()==b.Hash().
    */
    virtual bool IsSame(Node const& other) const;


    /**
    \brief Compute the derivative with respect to a single variable.

    If no variable is passed, produces a Jacobian tree when all is said and done, which is used for evaluating the Jacobian.

    Each distinct node is differentiated once per variable, however many paths reach it: the
    derivative of a subexpression is kept for as long as the outermost call (or an enclosing
    DerivativeScope) lasts, and every later request for it is answered from there.  A node type
    supplies its own rule as DifferentiateImpl, which calls this on its operands.

    \param v The variable to differentiate with respect to.
    \return The result of differentiation.  Jacobian or regular Node depending on what you passed in.
    */
    std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const;

    /**
    \brief Differentiate repeatedly with respect to a single variable.

    Applies the single-variable derivative `count` times, e.g. `Differentiate(x, 2)` is the
    second partial derivative with respect to `x`.  A `count` of zero returns the node itself.

    \param v The variable to differentiate with respect to.
    \param count How many times to differentiate.
    \return The `count`-th partial derivative with respect to `v` (a regular Node).
    */
    std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v, unsigned count) const;

    /**
    \brief Differentiate with respect to each variable in a group, in sequence.

    Applies the single-variable derivative once for each variable in `vars`, in order, e.g.
    `Differentiate({x, x, y})` is `∂³/∂y∂x²`.  Repeated entries are allowed.  An empty group
    returns the node itself.  (For polynomials mixed partials commute, so the order is immaterial.)

    \param vars The variables to differentiate with respect to, in application order.
    \return The mixed partial derivative (a regular Node).
    */
    std::shared_ptr<Node> Differentiate(VariableGroup const& vars) const;

    /**
    Compute the degree, optionally with respect to a single variable.

    Answered from the node's degree record (see DegreeRecord), which is made once per node from
    its operands' records, so asking costs a lookup however large the expansion beneath the node.
    The per-type rule lives in DegreeImpl.

    \param v Shared pointer to variable with respect to which you want to compute the degree of the Node.
    \return The degree.  Will be negative if the Node is non-polynomial.
    */
    int Degree(std::shared_ptr<Variable> const& v = nullptr) const;

    /// \brief Per-node-type degree in one variable.  Asks its operands through Degree(), so
    /// their records answer; a node's own record is made by calling this, so it must never ask
    /// Degree() of its own node.
    /// \param v The variable, or null for the total degree.
    /// \return The degree.  Negative if the node is non-polynomial.
    virtual int DegreeImpl(std::shared_ptr<Variable> const& v) const = 0;

    /**
     Compute the multidegree with respect to a variable group.

     Answered from the node's degree record, like Degree.  The per-type rule lives in
     MultiDegreeImpl.

     Why it matters (b2#417, #505): the node graph is hash-consed, so a small DAG can stand for
     an enormous expansion -- measured, a 46-node expression whose printed form is 196 KB.
     `CanonicalizeNaryOperands` asks this of every operand of every sum and product it builds.
     A plain recursive walk visits the TREE, so it cost O(expansion), which is how a cofactor
     expansion inside a deflation stage came to run for hours with a DAG that never exceeded
     1102 nodes.  A walk memoized within one call (#417) costs O(DAG), and asked once per
     construction that is still O(DAG^2) for building an expression: the derivative of a 7 x 7
     determinant took seconds.  The record makes each construction's question O(operands).

     #417 measured a cache of RESULTS between calls 29% slower than the memoized walk, because
     filling it cost the same walk.  The record does not walk: it is made from the operands'
     records, which already exist for any operand that has been built into a sum or product.
    */
    std::vector<int> MultiDegree(VariableGroup const& vars) const;

    /// \brief Per-node-type multidegree.  Asks its operands through MultiDegree(), and its own
    /// node only through DegreeImpl, never Degree() (see DegreeImpl).
    /// \param vars The variable group.
    /// \return One degree per variable of the group.
    virtual std::vector<int> MultiDegreeImpl(VariableGroup const& vars) const = 0;

    /**
    \brief The variables in this expression, distinct and sorted by name (GatherVariables).

    Part of the node's degree record: made once per node, from its operands' variables.

    \return The variables.
    */
    VariableGroup const& Variables() const;

    /**
    \brief This node's content key, the canonical encoding (CanonicalEncoding) of the tree it roots.

    Canonicalization orders the operands of a sum or product by it when their multidegrees tie.
    Computed on the first request and kept: a node is immutable, so its encoding never changes.

    \return The canonical encoding.
    */
    std::string const& ContentKey() const;

    /**
     Compute the overall degree with respect to a variable group.

    \param vars A group of variables.
     \return The degree.  Will be negative if the Node is non-polynomial.
    */
    virtual int Degree(VariableGroup const& vars) const = 0;

    /**
    Homogenize a tree, returning a NEW homogenized tree (functional / non-mutating).  Input a
    variable group holding the non-homogeneous variables, and the new homogenizing variable.
    The homvar may be an element of the variable group, that's perfectly ok.

    The input tree is never modified; degree-deficient summands are padded with powers of homvar
    in a freshly-built tree (so a throw on a non-polynomial term can't leave a half-homogenized
    tree behind).  The default (leaves, and anything with nothing to homogenize) returns the node
    unchanged.

    \param homvar The homogenizing variable, which is multiplied against terms with degree deficiency with repect to other terms.
    \param vars A group of variables, with respect to which you wish to homogenize.
    \return A homogenized tree (possibly the same node, if nothing changed).
    */
    virtual std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const;

    /**
    Check for homogeneity, absolutely with respect to all variables, including path variables and all other variable types, or with respect to a single varaible, if passed.

    \return True if it is homogeneous, false if not.
    */
    virtual bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const = 0;

    /**
    Check for homogeneity, with respect to a variable group.

    \return True if it is homogeneous, false if not.
    */
    virtual bool IsHomogeneous(VariableGroup const& vars) const = 0;

    ///////// PUBLIC PURE METHODS /////////////////

    /**
    Check if a Node is polynomial -- it has degree at least 0.  Negative degrees indicate non-polynomial.

    \return True if it is polynomial, false if not.
    */
    bool IsPolynomial(std::shared_ptr<Variable> const&v = nullptr) const;


    /**
    Check if a Node is polynomial -- it has degree at least 0.  Negative degrees indicate non-polynomial.

    \return True if it is polynomial, false if not.
    */
    bool IsPolynomial(VariableGroup const&v) const;






protected:
    /// Memoized structural hash (computed on first Hash() call; nodes are immutable so it
    /// never needs invalidating).  Not serialized -- a clone recomputes it on demand.
    mutable std::optional<std::size_t> structural_hash_;

    /// Opaque compiled-evaluator cache (a system-layer StraightLineProgram for [this]); see
    /// EvalProgram().  Lazily populated by EvalExpression; transient, not serialized.
    mutable std::shared_ptr<const void> eval_program_;

    /**
    \brief What Variables, Degree and MultiDegree report about this node, made once per node.

    Its variables; for each variable it mentions (its variables and its differentials'
    variables), its degree and multidegree entry in that variable; its total degree; and one
    degree and one multidegree entry for every variable it does not mention.  The last are a
    single pair because the degree rules compare a variable only with the variables a subtree
    mentions, so every unmentioned variable is answered alike -- not always 0: a negative
    constant power reports -1 in every variable.  Made from the operands' records, so it costs
    O(operands x mentioned variables), never a walk.  Defined in node.cpp.
    */
    struct DegreeRecord;

    /// The degree record, once made (see Record()).  Published once, atomically, and never
    /// replaced, so a reference into it lives as long as the node.  Transient, not serialized.
    mutable std::shared_ptr<const DegreeRecord> degree_record_;

    /// The content key, once computed (see ContentKey()).  Published like degree_record_.
    mutable std::shared_ptr<const std::string> content_key_;

    /// \brief This node's degree record, made on first request from its operands' records.
    /// \return The record.
    DegreeRecord const& Record() const;

    /// \brief A variable that is in no expression, for the record's answer about unmentioned variables.
    /// \return The variable.
    static std::shared_ptr<Variable> const& Unmentioned();

    /// Compute this node's structural hash.  Default: identity (the object address), so
    /// distinct nodes hash distinctly.  Value/operator nodes override to be structural.
    virtual std::size_t HashImpl() const;

    /**
    \brief This node type's differentiation rule, applied to this node alone.

    Called by Differentiate only when the derivative is not already known.  Differentiate an
    operand by calling its Differentiate, never its DifferentiateImpl, so that a subexpression
    shared by many operands is differentiated once.

    \param v The variable to differentiate with respect to.
    \return The derivative of this node.
    */
    virtual std::shared_ptr<Node> DifferentiateImpl(std::shared_ptr<Variable> const& v) const = 0;

    /// Combine an extra value into a running hash (boost::hash_combine recipe).
    static void HashCombine(std::size_t& seed, std::size_t value)
    {
        seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    }



    Node();

private:
    friend std::ostream& operator<<(std::ostream & out, const Node& N);

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, const unsigned /*version*/) {
        register_derived_node_types(ar);
    }

}; // re: class node

    /**
    \brief While one of these is alive, every derivative computed on this thread is kept and reused.

    Node::Differentiate opens one for the length of each outermost call, so a subexpression
    reached by many paths is differentiated once, not once per path: the sub-minors of a
    determinant built by cofactor expansion used to be differentiated k! times at the top (#505).
    Open one around a loop that differentiates many expressions, such as a whole Jacobian, so that
    subexpressions they share are differentiated once among all of them.  Scopes nest: an inner one
    uses the outer one's derivatives.  It keeps every node it has differentiated, and every
    derivative, alive until it closes.
    */
    class DerivativeScope
    {
    public:
        /// \brief Start keeping derivatives on this thread, unless an enclosing scope already is.
        DerivativeScope();
        /// \brief Stop keeping them and release them, if this scope started it.
        ~DerivativeScope();
        DerivativeScope(DerivativeScope const&) = delete;
        DerivativeScope& operator=(DerivativeScope const&) = delete;

        /// \cond DERIVATIVE_SCOPE_DETAIL
        struct Derivatives;
        /// \endcond
    private:
        std::unique_ptr<Derivatives> owned_;   ///< the derivatives, when this scope is the outermost
    };

    /**
     a single layer of indirection, though which to call the overridden virtual print() method which must be defined for each non-abstract Node type.
     */
    inline std::ostream& operator<<(std::ostream & out, const Node& N)
    {
        N.print(out);
        return out;
    }

    /// \brief Stream insertion for a Node pointer: dispatches to the node's virtual print().
    inline std::ostream& operator<<(std::ostream & out, const std::shared_ptr<Node>& N)
    {
        N->print(out);
        return out;
    }

    /**
    \brief Print a tree in one dialect (see print_dialect.hpp).

    \param N The tree to print.
    \param dialect Who the text is for.
    \return The printed expression.
    */
    inline std::string PrintIn(Node const& N, PrintDialect dialect)
    {
        std::ostringstream out;
        SetDialect(out, dialect);
        N.print(out);
        return out.str();
    }

    /// \brief The Bertini 1 spelling of a tree: `^` for powers, `(re+im*I)` for complex
    ///        constants, every constant with all its digits.  What the classic parser reads back.
    inline std::string PrintClassic(Node const& N) { return PrintIn(N, PrintDialect::Classic); }

    /**
    \brief The Python spelling of a tree: `**` for powers.

    \param N The tree to print.
    \param exact False for `str` (readable; the same constant shapes as Classic), true for
           `repr` (constant spellings `eval` rebuilds at full precision in the `bertini` namespace).
    \return The printed expression.
    */
    inline std::string PrintPython(Node const& N, bool exact)
    {
        return PrintIn(N, exact ? PrintDialect::PythonExact : PrintDialect::PythonReadable);
    }


    /**
    \brief Hash-cons a freshly-built node: return an existing structurally-equal node if one is
    live, otherwise register and return this one.

    The intern table is a process-global, weak (self-cleaning) map keyed by Node::Hash() and
    disambiguated by Node::IsSame().  Every Make() routes its just-constructed node through here,
    so structurally-equal subtrees collapse to a single shared object (hash-consing).  On a hit
    the just-built candidate is discarded.  Nodes whose IsSame() is identity (e.g. Variable,
    Function, Pi, E) never match, so they pass through unchanged -- no special-casing.

    Thread note: guarded by a mutex, contended only during single-threaded authoring.
    Deserialization constructs nodes WITHOUT going through Make/Intern; whatever holds the
    loaded nodes re-interns them as the last step of its load (System and Slice do, ADR-0068),
    so nothing a caller receives is outside this table.
    */
    std::shared_ptr<Node> Intern(std::shared_ptr<Node> const& candidate);


    } // re: namespace node
} // re: namespace bertini



#endif
/* defined(BERTINI_NODE_BASE_HPP) */
