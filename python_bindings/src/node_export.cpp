//This file is part of Bertini 2.
//
//python/node_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/node_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/node_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  James Collins
//  West Texas A&M University
//  Spring 2016
//
//  silviana amethyst
//  UWEC
//  2017, Spring 2018
//
//
//  python/node_export.cpp:  Source file for exposing Node class to python.

#include <stdio.h>
#include <sstream>


#include <boost/python/raw_function.hpp>

#include "node_export.hpp"
#include "bertini2/function_tree/canonical.hpp"
#include "bertini2/system/eval_expression.hpp"


namespace bertini{
	namespace python{
		
		// Wrapper struct to allow derived classes to overide methods in python
		struct NodeWrap : Node, wrapper<Node>
		{
			int Degree(std::shared_ptr<Variable> const& v = nullptr) const {return this->get_override("Degree")(v); }
			int Degree(VariableGroup const& vars) const {return this->get_override("Degree")(vars); }
			
			std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const {return this->get_override("Differentiate")(v); }
			
			std::vector<int> MultiDegree(VariableGroup const& vars) const {return this->get_override("MultiDegree")(vars); }
			
			std::shared_ptr<Node> Homogenized(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) const { return this->get_override("Homogenized")(vars, homvar); }
			
			bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const {return this->get_override("IsHomogeneous")(v); }
			bool IsHomogeneous(VariableGroup const& vars) const {return this->get_override("IsHomogeneous")(vars); }
			
			bool IsPolynomial(std::shared_ptr<Variable> const&v = nullptr) const {return this->get_override("IsPolynomial")(v); }
			bool IsPolynomial(VariableGroup const&v) const {return this->get_override("IsPolynomial")(v); }
			

			
		}; // re: NodeWrap





		
	
		
		
		
		// Python's power operator is ** (^ is bitwise-xor), but the function tree prints powers with
		// ^.  __repr__ should be copy-pasteable Python, so translate ^ -> ** (^ only ever denotes a
		// power in the tree's printed form).  __str__ keeps the classic ^ (which the input parser reads).
		static std::string NodeRepr(std::shared_ptr<Node> const& self)
		{
			std::ostringstream oss;
			oss << *self;
			std::string const s = oss.str();
			std::string out;
			out.reserve(s.size() + 8);
			for (char c : s)
			{
				if (c == '^') out += "**";
				else          out += c;
			}
			return out;
		}

		template<typename NodeBaseT>
		template<class PyClass>
		void NodeVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("degree", &Deg0, (arg("self")),"compute the algebraic degree of node in a function tree, with respect to all variables. returns one integer.  negative is non-algebraic.")
			.def("degree", Deg1, (arg("self"),arg("var")),"compute the algebraic degree of node in a function tree, with respect to a particular variable. returns one integer.  negative is non-algebraic.")
			.def("degree", Deg2, (arg("self"),arg("vars")),"compute the algebraic degree of node in a function tree, with respect to a variable group. returns one integer.  negative is non-algebraic.")
			.def("differentiate", Diff0, (arg("self")),"differentiate a node.  is with respect to all variables.  you get a Jacobian back, which represents derivatives wrt all variables simultaneously.")
			.def("differentiate", Diff1, (arg("self")),"differentiate a node with respect to one variable.  You get a regular old Node in a Function Tree back.")
			.def("differentiate", DiffN, (arg("self"),arg("var"),arg("count")),"differentiate a node with respect to one variable, `count` times (e.g. n.differentiate(x, 2) for the second partial).  count=0 returns the node itself.  You get a regular Node back.")
			.def("differentiate", DiffList, (arg("self"),arg("vars")),"differentiate a node with respect to each variable in `vars`, in sequence -- a mixed partial (e.g. n.differentiate([x, x, y])).  Repeated entries are allowed; an empty list returns the node itself.  You get a regular Node back.")
			.def("multidegree", &NodeBaseT::MultiDegree, (arg("self"),arg("vars")),"Compute an integer vector containing the degrees with respect to the variables in `vars`.  Negative entries indicate non-polynomiality")
			.def("homogenized", &NodeBaseT::Homogenized, (arg("self"),arg("vars"), arg("homvar")), "Return a NEW homogenized copy of this function tree (non-mutating) with respect to the variables in `vars` using the homogenizing variable `homvar`.  Degree-deficient terms are padded with powers of `homvar` so all terms share the same degree.  The original tree is left untouched.")
			.def("is_homogeneous", IsHom0,(arg("self")), "test if this Node is homogeneous with respect to all Variables.")
			.def("is_homogeneous", IsHom1,(arg("self"),arg("var")), "test if this Node is homogeneous with respect to Variable `var`.")
			.def("is_homogeneous", IsHom2,(arg("self"),arg("vars")), "test if this Node is homogeneous with respect to the Variables in `vars`.")
			.def("is_polynomial", IsPoly0,(arg("self")), "test if this Node is polynomial with respect to all Variables.")
			.def("is_polynomial", IsPoly1,(arg("self"),arg("var")), "test if this Node is polynomial with respect to Variable `var`.")
			.def("is_polynomial", IsPoly2,(arg("self"),arg("vars")), "test if this Node is polynomial with respect to Variables `vars`.")

			.def(self_ns::str(self_ns::self))
			.def("__repr__", &NodeRepr)
			
			.def("__add__",addNodeNode)
			.def("__add__",addNodeMpfr)
			.def("__radd__",&raddNodeMpfr)
			.def("__add__",addNodeRat)
			.def("__radd__",&raddNodeRat)
			.def("__add__",addNodeInt)
			.def("__radd__",raddNodeInt)
			.def("__iadd__",&NodeVisitor::iaddNodeNode)
			.def("__iadd__", &NodeVisitor::iaddSumNode)
			
			.def("__sub__",subNodeNode)
			.def("__sub__",subNodeMpfr)
			.def("__rsub__",rsubNodeMpfr)
			.def("__sub__",subNodeRat)
			.def("__rsub__",rsubNodeRat)
			.def("__sub__",subNodeInt)
			.def("__rsub__",rsubNodeInt)
			.def("__isub__",&NodeVisitor::isubNodeNode)
			.def("__isub__", &NodeVisitor::isubSumNode)
			
			.def("__mul__",multNodeNode)
			.def("__mul__",multNodeMpfr)
			.def("__rmul__",rmultNodeMpfr)
			.def("__mul__",multNodeRat)
			.def("__rmul__",rmultNodeRat)
			.def("__mul__",multNodeInt)
			.def("__rmul__",rmultNodeInt)
			.def("__imul__",&NodeVisitor::imultNodeNode)
			.def("__imul__",imultMultNode)
			

			.def("__div__",divNodeNode)
			.def("__truediv__",divNodeNode)
			.def("__itruediv__",&NodeVisitor::idivNodeNode)



			
			.def("__div__",divNodeMpfr)
			.def("__truediv__",divNodeMpfr)

			.def("__rdiv__",rdivNodeMpfr)
			.def("__rtruediv__",rdivNodeMpfr)

			.def("__rdiv__",rdivNodeRat)
			.def("__rtruediv__",rdivNodeRat)

			.def("__div__",divNodeInt)
			.def("__truediv__",divNodeInt)

			.def("__rdiv__",rdivNodeInt)
			.def("__rtruediv__",rdivNodeInt)

			.def("__idiv__",&NodeVisitor::idivNodeNode)
			.def("__itruediv__",&NodeVisitor::idivNodeNode)
			
			.def("__idiv__",idivMultNode)
			.def("__itruediv__",idivMultNode)
			
			.def("__neg__", negNode)
			
			.def("__pow__",powNodeNode)
			.def("__pow__",powNodeMpfr)
			.def("__pow__",powNodeRat)
			.def("__pow__",powNodeInt)
			;
			
			
			def("exp", expNodeNode, "the symbolic exponential operator");
			def("log", logNodeNode, "the symbolic natural log operator");
			def("sin", sinNodeNode, "the symbolic sine operator");
			def("asin", asinNodeNode, "the symbolic arcsine operator");
			def("cos", cosNodeNode, "the symbolic cosine operator");
			def("acos", acosNodeNode, "the symbolic arccosine operator");
			def("tan", tanNodeNode, "the symbolic tangent operator");
			def("atan", atanNodeNode, "the symbolic arctangent operator");

			// ---- canonical operand ordering ----
			enum_<node::MonomialOrder>("MonomialOrder",
				"the monomial order used to canonically order Sum/Mult operands")
				.value("Lex",     node::MonomialOrder::Lex)
				.value("RevLex",  node::MonomialOrder::RevLex)
				.value("GrevLex", node::MonomialOrder::GrevLex)
				;

			bool (*canon_get)()     = &node::CanonicalizeByDefault;
			void (*canon_set)(bool) = &node::SetCanonicalizeByDefault;
			def("canonicalize", canon_get,
				"whether Sum/Mult operands are canonically ordered, so x+y and y+x are one node");
			def("canonicalize", canon_set, (arg("on")),
				"enable/disable canonical operand ordering, session-global (the per-expression opt-out)");

			node::MonomialOrder (*order_get)()                   = &node::CurrentMonomialOrder;
			void                (*order_set)(node::MonomialOrder) = &node::SetMonomialOrder;
			def("monomial_order", order_get, "the current monomial order used for canonicalization");
			def("monomial_order", order_set, (arg("order")),
				"set the monomial order (Lex/RevLex/GrevLex), session-global");

		}

		
		
		
		// Interpret a single Python value (int, float, complex, or a multiprecision
		// number) as an complex_mp.  Multiprecision inputs are taken as-is; native
		// Python numbers are embedded at the current default precision (so a Python
		// float carries only float64 worth of information --- the documented cap).
		static complex_mp CoerceToMpfrComplex(object const& o)
		{
			extract<complex_mp> as_mp(o);
			if (as_mp.check()) return as_mp();

			extract<std::complex<double>> as_complex(o);
			if (as_complex.check()) { auto c = as_complex(); return complex_mp(c.real(), c.imag()); }

			extract<double> as_double(o);
			if (as_double.check()) return complex_mp(as_double());

			extract<long> as_long(o);
			if (as_long.check()) return complex_mp(static_cast<double>(as_long()));

			throw std::runtime_error("could not interpret a supplied value as a number in eval");
		}

		// Resolve a Python object naming a variable -- a Variable node or a name string -- to its name.
		static std::string VariableNameOf(object const& key)
		{
			extract<std::string> as_str(key);
			if (as_str.check())
				return as_str();
			extract<std::shared_ptr<node::Variable>> as_var(key);
			if (as_var.check())
				return as_var()->name();
			throw std::runtime_error("eval: dictionary keys must be Variables or variable-name strings");
		}

		// Parse a Python object naming a variable ordering -- a VariableGroup or a sequence of
		// Variable nodes -- into a VariableGroup.
		static VariableGroup OrderingOf(object const& vobj)
		{
			extract<VariableGroup> as_vg(vobj);
			if (as_vg.check())
				return as_vg();
			VariableGroup vars;
			for (long i = 0; i < len(vobj); ++i)
				vars.push_back(extract<std::shared_ptr<node::Variable>>(vobj[i])());
			return vars;
		}

		// f.eval(...) --- evaluate this expression at a point.  No System is required; values bind by
		// variable name and the result is a complex_mp.  Call forms (issue #300):
		//   (a) keyword arguments naming the variables:  f.eval(x=2, y=5)
		//   (b) a single positional dict {Variable-or-name: value}:  f.eval({x: 2, y: 5})
		//   (c) a 1-D array / list, mapped in order to the expression's variables() (which are sorted
		//       by name); override that ordering with a second positional argument or a variables=
		//       keyword:  f.eval(pt)  /  f.eval(pt, [x, y, z])  /  f.eval(pt, variables=[x, y, z])
		// A missing value for a variable of the expression always throws.  By default (strict=False)
		// supplied names that are not variables of the expression are ignored (the expression is
		// constant with respect to them) -- so one full point/ordering can be reused across an
		// expression and its derivatives; pass strict=True to reject such extras as likely typos.
		static object NodeEvalRaw(tuple args, dict kwargs)
		{
			long const nargs = len(args);
			if (nargs < 1)
				throw std::runtime_error("eval: missing self");
			std::shared_ptr<Node> self = extract<std::shared_ptr<Node>>(args[0]);

			bool strict = false;   // default: ignore variables the expression does not depend on
			if (kwargs.has_key("strict"))
				strict = extract<bool>(kwargs["strict"]);

			std::map<std::string, complex_mp> values;

			// (a) keyword form -- f.eval(x=2, y=5)
			if (nargs == 1)
			{
				list items = dict(kwargs).items();
				for (long i = 0; i < len(items); ++i)
				{
					object pair = items[i];
					std::string name = extract<std::string>(pair[0]);
					if (name == "variables" || name == "strict")   // control kwargs, not variable values
						continue;
					values[name] = CoerceToMpfrComplex(object(pair[1]));
				}
				return object(bertini::EvalExpression<complex_mp>(self, values, strict));
			}

			if (nargs > 3)
				throw std::runtime_error("eval: pass a point (a dict, or a 1-D array/list) and at most "
				                         "an ordering; got too many positional arguments");

			object point = args[1];

			// (b) a dict {Variable-or-name: value}
			extract<dict> as_dict(point);
			if (as_dict.check())
			{
				if (nargs > 2)
					throw std::runtime_error("eval: a variable ordering is meaningful only for a 1-D "
					                         "array/list point, not a dict");
				dict d = as_dict();
				list items = d.items();
				for (long i = 0; i < len(items); ++i)
				{
					object pair = items[i];
					values[VariableNameOf(object(pair[0]))] = CoerceToMpfrComplex(object(pair[1]));
				}
				return object(bertini::EvalExpression<complex_mp>(self, values, strict));
			}

			// (c) a 1-D array/list mapped to a variable ordering.  The ordering may be given as a
			// second positional argument or a variables= keyword; otherwise it defaults to the
			// expression's own variables() (sorted by name).
			VariableGroup vars;
			if (nargs > 2)
				vars = OrderingOf(args[2]);
			else if (kwargs.has_key("variables"))
				vars = OrderingOf(kwargs["variables"]);
			else
				vars = bertini::node::GatherVariables(self);   // the expression's own variables, sorted by name

			long const n = len(point);
			if (static_cast<size_t>(n) != vars.size())
				throw std::runtime_error("eval: the point has " + std::to_string(n) + " entries but the "
				                         "ordering has " + std::to_string(vars.size()) + " variables; "
				                         "pass a matching ordering if they differ");
			for (long i = 0; i < n; ++i)
				values[vars[static_cast<size_t>(i)]->name()] = CoerceToMpfrComplex(object(point[i]));

			return object(bertini::EvalExpression<complex_mp>(self, values, strict));
		}

		// f.variables() --- the distinct variables appearing in this expression, sorted by name.
		static VariableGroup NodeVariables(std::shared_ptr<Node> const& self)
		{
			return bertini::node::GatherVariables(self);
		}

		// f.simplify() --- a NEW, functionally simplified copy (non-mutating): literal zeros/ones
		// vanish, exact constants fold (including constant powers, so (3**2)->9 and i**2->-1).
		static std::shared_ptr<Node> NodeSimplify(std::shared_ptr<Node> const& self)
		{
			return self->Simplified();
		}

		// Coerce a Python object used as a substitution value into a Node.  A Node (a Variable, or
		// any expression) is used as-is; an exact Python int becomes an Integer; a fractions.Fraction
		// becomes an exact Rational; anything else numeric becomes a Complex (mpfr, possibly lossy).
		static std::shared_ptr<Node> CoerceSubValueToNode(object const& o)
		{
			extract<std::shared_ptr<Node>> as_node(o);
			if (as_node.check()) return as_node();

			if (PyLong_Check(o.ptr()))
			{
				extract<int> as_int(o);
				if (as_int.check())
					return node::Integer::Make(as_int());
				return node::Integer::Make(std::string(extract<std::string>(str(o))));  // big int, via decimal string
			}

			object Fraction = import("fractions").attr("Fraction");
			if (PyObject_IsInstance(o.ptr(), Fraction.ptr()) == 1)
			{
				int num = extract<int>(o.attr("numerator"));
				int den = extract<int>(o.attr("denominator"));
				return node::Rational::Make(num, den, 0, 1);   // real rational num/den, zero imaginary part
			}

			return node::Complex::Make(CoerceToMpfrComplex(o));
		}

		// f.subs(...) --- symbolically substitute variables, returning a NEW expression (a Node).  The
		// substitution is simultaneous and single-pass; a variable not named is left unchanged.  Forms:
		//   (a) keyword arguments:      f.subs(x=3, y=5)
		//   (b) a dict {var-or-name: value}:  f.subs({x: 3, y: z**2})
		//   (c) a variable and its value:     f.subs(x, z**2)
		// Values may be Nodes (a Variable or any expression) or numbers (int -> Integer, Fraction ->
		// Rational, else Complex).
		static object NodeSubsRaw(tuple args, dict kwargs)
		{
			long const nargs = len(args);
			if (nargs < 1)
				throw std::runtime_error("subs: missing self");
			std::shared_ptr<Node> self = extract<std::shared_ptr<Node>>(args[0]);

			SubstitutionMap subs;

			if (nargs == 1)   // (a) keyword form
			{
				list items = dict(kwargs).items();
				for (long i = 0; i < len(items); ++i)
				{
					object pair = items[i];
					subs[extract<std::string>(pair[0])] = CoerceSubValueToNode(object(pair[1]));
				}
			}
			else if (nargs == 2)   // (b) a dict
			{
				extract<dict> as_dict(args[1]);
				if (!as_dict.check())
					throw std::runtime_error("subs: a single positional argument must be a dict "
					                         "{variable: value}; for one substitution use subs(var, value)");
				dict d = as_dict();
				list items = d.items();
				for (long i = 0; i < len(items); ++i)
				{
					object pair = items[i];
					subs[VariableNameOf(object(pair[0]))] = CoerceSubValueToNode(object(pair[1]));
				}
			}
			else if (nargs == 3)   // (c) variable, value
			{
				subs[VariableNameOf(args[1])] = CoerceSubValueToNode(object(args[2]));
			}
			else
				throw std::runtime_error("subs: pass keyword arguments, a dict {variable: value}, "
				                         "or a variable and its value");

			return object(self->Subs(subs));
		}

		void ExportNode()
		{
			class_<NodeWrap, boost::noncopyable, Nodeptr >("AbstractNode", no_init)
			.def(NodeVisitor<Node>())
			.def("eval", raw_function(&NodeEvalRaw, 1),
				"evaluate this expression at a point, returning a complex_mp.  No System is needed.  "
				"Forms: keyword args f.eval(x=2, y=5); a dict f.eval({x: 2, y: 5}) (keys may be Variables "
				"or name strings); or a 1-D array/list f.eval(pt) mapped in order to the expression's "
				"variables() (sorted by name) -- override the ordering with a second positional argument "
				"or variables=, e.g. f.eval(pt, [x, y, z]) or f.eval(pt, variables=[x, y, z]).  "
				"By default (strict=False) supplied variables the expression does not depend on are "
				"ignored (it is constant with respect to them), so one full point can be reused across "
				"an expression and its derivatives; pass strict=True to reject such extras as typos.  A "
				"missing value for a variable the expression DOES depend on always raises.  Evaluation is "
				"at the current default precision; native Python floats carry only float64 of information.")
			.def("variables", &NodeVariables, (arg("self")),
				"The distinct variables appearing in this expression, sorted by name.")
			.def("simplify", &NodeSimplify, (arg("self")),
				"Return a NEW, functionally simplified copy of this expression (non-mutating): literal "
				"zeros/ones vanish and exact constants fold, including constant powers -- so (3**2) "
				"simplifies to 9 and i**2 (i = Complex(0,1)) to -1.  The original is left untouched.")
			.def("subs", raw_function(&NodeSubsRaw, 1),
				"symbolically substitute variables, returning a NEW expression.  Forms: keyword args "
				"f.subs(x=3, y=5); a dict f.subs({x: 3, y: z**2}) (keys may be Variables or name strings); "
				"or a single variable and its value f.subs(x, z**2).  Values may be nodes (a Variable or "
				"any expression) or numbers (int -> Integer, fractions.Fraction -> Rational, else Complex).  "
				"Substitution is simultaneous and single-pass -- {x: y, y: z} maps x+y to y+z (no cascade), "
				"and the order is irrelevant; a variable not named is left unchanged.  Only variables are "
				"substituted (not subexpressions), and the result is simplified.")
			;
		};
		
		
	} //namespace python
} // namespace bertini
