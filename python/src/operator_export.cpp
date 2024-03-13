//This file is part of Bertini 2.
//
//python/operator_export.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/operator_export.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/operator_export.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2016-2018 by Bertini2 Development Team
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
//  Spring 2018
//
//
//
//  python/operator_export.cpp:  Source file for exposing operator nodes to python.



#include <stdio.h>
#include "operator_export.hpp"



namespace bertini{
	namespace python{

		struct UnaryOpWrap : UnaryOperator, wrapper<UnaryOperator>
		{
			void SetOperand(std::shared_ptr<Node> new_child)
			{
				if (override SetOperand = this->get_override("SetOperand"))
					SetOperand(new_child); // *note*
				
				UnaryOperator::SetOperand(new_child);
			}
			void default_SetChild(std::shared_ptr<Node> new_child){ return this->UnaryOperator::SetOperand(new_child);}
		}; // re: NodeWrap

		
		struct NaryOpWrap : NaryOperator, wrapper<NaryOperator>
		{
			void AddOperand(std::shared_ptr<Node> child)
			{
				if (override AddOperand = this->get_override("AddOperand"))
					AddOperand(child); // *note*
				
				NaryOperator::AddOperand(child);
			}
			void default_AddOperand(std::shared_ptr<Node> child){ return this->NaryOperator::AddOperand(child);}
		}; // re: NodeWrap

		
		
		
		template<typename NodeBaseT>
		template<class PyClass>
		void UnaryOpVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("set_operand", &NodeBaseT::SetOperand )
			.def("operand", &NodeBaseT::Operand )
			;
		}
		
		template<typename NodeBaseT>
		template<class PyClass>
		void NaryOpVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("add_operand", &NodeBaseT::AddOperand )
			.def("first_operand", &NodeBaseT::FirstOperand )
			.def("num_operands", &NodeBaseT::NumOperands )
			;
		}

		template<typename NodeBaseT>
		template<class PyClass>
		void SumMultOpVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			
			cl
			.def("add_operand", AddOperand2)
			;
		}


		template<typename NodeBaseT>
		template<class PyClass>
		void PowerOpVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("set_exponent", &PowerOperator::SetExponent)
			.def("set_base", &PowerOperator::SetBase)
			;
		}
		
		template<typename NodeBaseT>
		template<class PyClass>
		void IntPowOpVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.add_property("exponent", getexp, setexp)
			;
		}



		void ExportOperators()
		{
			scope current_scope;
			std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
			new_submodule_name.append(".operator");
			object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
			current_scope.attr("operator") = new_submodule;

			scope new_submodule_scope = new_submodule;

			// Operator class
			class_<Operator, boost::noncopyable, bases<Node>, std::shared_ptr<Operator> >("AbstractOp", no_init)
			;

			// UnaryOperator class
			class_<UnaryOpWrap, boost::noncopyable, bases<Operator>, std::shared_ptr<UnaryOperator> >("Unary", no_init)
			.def(UnaryOpVisitor<UnaryOperator>())
			;

			// NaryOperator class
			class_<NaryOpWrap, boost::noncopyable, bases<Operator>, std::shared_ptr<NaryOperator> >("Nary", no_init)
			.def(NaryOpVisitor<NaryOperator>())
			;

			// SumOperator class
			class_<SumOperator, bases<NaryOperator>, std::shared_ptr<SumOperator> >("Sum", no_init)
			.def("__init__", make_constructor(&SumOperator::template Make<const Nodeptr&, const Nodeptr &> ))
			.def("__init__", make_constructor(&SumOperator::template Make<const Nodeptr&, bool const&, const Nodeptr&, bool const&> ))
			.def(SumMultOpVisitor<SumOperator>())
			;


			// NegateOperator class
			class_<NegateOperator, bases<UnaryOperator>, std::shared_ptr<NegateOperator> >("Negate", no_init )
			.def("__init__", make_constructor(&NegateOperator::template Make<const Nodeptr&>))
			;

			// MultOperator class
			class_<MultOperator, bases<NaryOperator>, std::shared_ptr<MultOperator> >("Mult", no_init )
			.def("__init__", make_constructor(&MultOperator::template Make<const Nodeptr&, const Nodeptr &>))
			.def("__init__", make_constructor(&MultOperator::template Make<const Nodeptr&, bool const&, const Nodeptr&, bool const&>))
			.def(SumMultOpVisitor<MultOperator>())
			;
			
			// PowerOperator class
			class_<PowerOperator, bases<Operator>, std::shared_ptr<PowerOperator> >("Power", no_init )
			.def("__init__", make_constructor(&PowerOperator::template Make<const Nodeptr&, const Nodeptr &>))
			.def(PowerOpVisitor<PowerOperator>())
			;
			
			// IntegerPowerOperator class
			class_<IntegerPowerOperator, bases<UnaryOperator>, std::shared_ptr<IntegerPowerOperator> >("IntegerPower", no_init )
			.def("__init__", make_constructor(&IntegerPowerOperator::template Make<const Nodeptr&, int const&>))
			.def(PowerOpVisitor<IntegerPowerOperator>())
			;

			// SqrtOperator class
			class_<SqrtOperator, bases<UnaryOperator>, std::shared_ptr<SqrtOperator> >("Sqrt", no_init)
			.def("__init__", make_constructor(&SqrtOperator::template Make<const Nodeptr&>))
			;

			// ExpOperator class
			class_<ExpOperator, bases<UnaryOperator>, std::shared_ptr<ExpOperator> >("Exp", no_init)
			.def("__init__", make_constructor(&ExpOperator::template Make<const Nodeptr&> ))
			;

			// LogOperator class
			class_<LogOperator, bases<UnaryOperator>, std::shared_ptr<LogOperator> >("Log", no_init)
			.def("__init__", make_constructor(&LogOperator::template Make<const Nodeptr&> ))
			;
			
			

			// TrigOperator class
			class_<TrigOperator, boost::noncopyable, bases<Node>, std::shared_ptr<TrigOperator> >("Trig", no_init)
			;

			// SinOperator class
			class_<SinOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<SinOperator> >("Sin", no_init)
			.def("__init__", make_constructor(&SinOperator::template Make<const Nodeptr&>))
			;
			// CosOperator class
			class_<CosOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<CosOperator> >("Cos", no_init)
			.def("__init__", make_constructor(&CosOperator::template Make<const Nodeptr&>))
			;
			// TanOperator class
			class_<TanOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<TanOperator> >("Tan", no_init)
			.def("__init__", make_constructor(&TanOperator::template Make<const Nodeptr&>))
			;

			// ArcSinOperator class
			class_<ArcSinOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<ArcSinOperator> >("ArcSin", no_init)
			.def("__init__", make_constructor(&ArcSinOperator::template Make<const Nodeptr&>))
			;
			// ArcCosOperator class
			class_<ArcCosOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<ArcCosOperator> >("ArcCos", no_init)
			.def("__init__", make_constructor(&ArcCosOperator::template Make<const Nodeptr&>))
			;
			// ArcTanOperator class
			class_<ArcTanOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<ArcTanOperator> >("ArcTan", no_init)
			.def("__init__", make_constructor(&ArcTanOperator::template Make<const Nodeptr&>))
			;


			
			
		}
	} //namespace python
} // namespace bertini



