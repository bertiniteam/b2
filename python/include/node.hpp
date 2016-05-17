//This file is part of Bertini 2.
//
//python/node.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/node.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/node.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//
//  Jeb Collins
//  West Texas A&M University
//  Mathematics
//  Fall 2015
//
//
//  python/node.hpp:  the header file for the python interface for Node class.


#ifndef BERTINI_PYTHON_NODE_HPP
#define BERTINI_PYTHON_NODE_HPP

#include "function_tree.hpp"


namespace bertini{
	namespace python{
		
		using namespace boost::python;
		using namespace bertini::node;
		
		using Node = node::Node;
		using NodePtr = std::shared_ptr<node::Node>;
		
		struct NodeWrap : Node, wrapper<Node>
		{
			void Reset()
			{
				if (override Reset = this->get_override("Reset"))
					Reset(); // *note*
				
				Node::Reset();
			}
			
			void default_Reset()
			{
				return this->Node::Reset();
			}
			
			void print(std::ostream& target) const
			{
				this->get_override("print")();
			}
		}; // re: NodeWrap

		// Addition operators
		std::shared_ptr<node::Node>(*addNodeNode)(std::shared_ptr<node::Node>, const std::shared_ptr<node::Node>&) = &operator+;
		std::shared_ptr<node::Node>(*addNodeDouble)(std::shared_ptr<node::Node>, double) = &operator+;
		std::shared_ptr<node::Node>(*addNodeDbl)(std::shared_ptr<node::Node>, dbl) = &operator+;
		std::shared_ptr<node::Node>(*addNodeMpfr)(std::shared_ptr<node::Node>, mpfr) = &operator+;
		std::shared_ptr<node::Node>(*addNodeInt)(std::shared_ptr<node::Node>, int) = &operator+;
		inline std::shared_ptr<node::Node> iaddNodeNode(std::shared_ptr<node::Node>  lhs, const std::shared_ptr<node::Node> & rhs)
		{
			return lhs += rhs;
		}
		inline std::shared_ptr<node::Node> iaddNodeDouble(std::shared_ptr<node::Node>  lhs, double rhs)
		{
			return lhs += rhs;
		}
		inline node::SumOperator iaddSumNode(node::SumOperator  lhs, const std::shared_ptr<node::Node> & rhs)
		{
			return lhs += rhs;
		}
		
		// Subtraction operators
		std::shared_ptr<node::Node>(*subNodeNode)(std::shared_ptr<node::Node>, const std::shared_ptr<node::Node>&) = &operator-;
		std::shared_ptr<node::Node>(*subNodeDouble)(std::shared_ptr<node::Node>, double) = &operator-;
		std::shared_ptr<node::Node>(*subNodeDbl)(std::shared_ptr<node::Node>, dbl) = &operator-;
		std::shared_ptr<node::Node>(*subNodeMpfr)(std::shared_ptr<node::Node>, mpfr) = &operator-;
		std::shared_ptr<node::Node>(*subNodeInt)(std::shared_ptr<node::Node>, int) = &operator-;
		inline std::shared_ptr<node::Node> isubNodeNode(std::shared_ptr<node::Node>  lhs, const std::shared_ptr<node::Node> & rhs)
		{
			return lhs -= rhs;
		}
		inline std::shared_ptr<node::Node> isubNodeDouble(std::shared_ptr<node::Node>  lhs, double rhs)
		{
			return lhs -= rhs;
		}
		inline node::SumOperator isubSumNode(node::SumOperator  lhs, const std::shared_ptr<node::Node> & rhs)
		{
			return lhs -= rhs;
		}
		
		
		// Negate operator
		std::shared_ptr<node::Node>(*negNode)(const std::shared_ptr<Node> &) = &operator-;

		
		
		
		// Multiplication operators
		std::shared_ptr<node::Node>(*multNodeNode)(std::shared_ptr<node::Node>, const std::shared_ptr<node::Node>&) = &operator*;
		std::shared_ptr<node::Node>(*multNodeDouble)(std::shared_ptr<node::Node>, double) = &operator*;
		std::shared_ptr<node::Node>(*multNodeDbl)(std::shared_ptr<node::Node>, dbl) = &operator*;
		std::shared_ptr<node::Node>(*multNodeMpfr)(std::shared_ptr<node::Node>, mpfr) = &operator*;
		std::shared_ptr<node::Node>(*multNodeInt)(std::shared_ptr<node::Node>, int) = &operator*;
		inline std::shared_ptr<node::Node> imultNodeNode(std::shared_ptr<node::Node>  lhs, const std::shared_ptr<node::Node> & rhs)
		{
			return lhs *= rhs;
		}
		inline std::shared_ptr<node::Node> imultNodeDouble(std::shared_ptr<node::Node>  lhs, double rhs)
		{
			return lhs *= rhs;
		}
		std::shared_ptr<node::Node>(*imultMultNode)(std::shared_ptr<node::MultOperator> &, const std::shared_ptr<node::Node> &) = &operator*=;

		// Division operators
		std::shared_ptr<node::Node>(*divNodeNode)(std::shared_ptr<node::Node>, const std::shared_ptr<node::Node>&) = &operator/;
		std::shared_ptr<node::Node>(*divNodeDouble)(std::shared_ptr<node::Node>, double) = &operator/;
		std::shared_ptr<node::Node>(*divNodeDbl)(std::shared_ptr<node::Node>, dbl) = &operator/;
		std::shared_ptr<node::Node>(*divNodeMpfr)(std::shared_ptr<node::Node>, mpfr) = &operator/;
		std::shared_ptr<node::Node>(*divNodeInt)(std::shared_ptr<node::Node>, int) = &operator/;
		inline std::shared_ptr<node::Node> idivNodeNode(std::shared_ptr<node::Node>  lhs, const std::shared_ptr<node::Node> & rhs)
		{
			return lhs /= rhs;
		}
		inline std::shared_ptr<node::Node> idivNodeDouble(std::shared_ptr<node::Node>  lhs, double rhs)
		{
			return lhs /= rhs;
		}
		std::shared_ptr<node::Node>(*idivMultNode)(std::shared_ptr<node::MultOperator> &, const std::shared_ptr<node::Node> &) = &operator/=;

		// Power operators
		std::shared_ptr<node::Node>(*powNodeNode)(const std::shared_ptr<node::Node> &, const std::shared_ptr<node::Node>&) = &pow;
		std::shared_ptr<node::Node>(*powNodeDouble)(const std::shared_ptr<node::Node>&, double) = &pow;
		std::shared_ptr<node::Node>(*powNodeDbl)(const std::shared_ptr<node::Node>&, dbl) = &pow;
		std::shared_ptr<node::Node>(*powNodeMpfr)(const std::shared_ptr<node::Node>&, mpfr) = &pow;
		std::shared_ptr<node::Node>(*powNodeInt)( std::shared_ptr<node::Node> const&, int) = &pow;
		
		// Transcindental operators
		std::shared_ptr<node::Node>(*expNodeNode)(const std::shared_ptr<node::Node> &) = &exp;
		std::shared_ptr<node::Node>(*logNodeNode)(const std::shared_ptr<node::Node> &) = &log;
		std::shared_ptr<node::Node>(*sinNodeNode)(const std::shared_ptr<node::Node> &) = &sin;
		std::shared_ptr<node::Node>(*asinNodeNode)(const std::shared_ptr<node::Node> &) = &asin;
		std::shared_ptr<node::Node>(*cosNodeNode)(const std::shared_ptr<node::Node> &) = &cos;
		std::shared_ptr<node::Node>(*acosNodeNode)(const std::shared_ptr<node::Node> &) = &acos;
		std::shared_ptr<node::Node>(*tanNodeNode)(const std::shared_ptr<node::Node> &) = &tan;
		std::shared_ptr<node::Node>(*atanNodeNode)(const std::shared_ptr<node::Node> &) = &atan;

		
	} // re: python
} // re: bertini


#endif
