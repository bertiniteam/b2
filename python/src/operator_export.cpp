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
// Copyright(C) 2016 by Bertini2 Development Team
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
//
//
//  python/operator_export.cpp:  Source file for exposing operator nodes to python.



#include <stdio.h>
#include "operator_export.hpp"



namespace bertini{
	namespace python{

		struct UnaryOpWrap : UnaryOperator, wrapper<UnaryOperator>
		{
			void SetChild(std::shared_ptr<Node> new_child)
			{
				if (override SetChild = this->get_override("SetChild"))
					SetChild(new_child); // *note*
				
				UnaryOperator::SetChild(new_child);
			}
			void default_SetChild(std::shared_ptr<Node> new_child){ return this->UnaryOperator::SetChild(new_child);}
		}; // re: NodeWrap

		
		struct NaryOpWrap : NaryOperator, wrapper<NaryOperator>
		{
			void AddChild(std::shared_ptr<Node> child)
			{
				if (override AddChild = this->get_override("AddChild"))
					AddChild(child); // *note*
				
				NaryOperator::AddChild(child);
			}
			void default_AddChild(std::shared_ptr<Node> child){ return this->NaryOperator::AddChild(child);}
		}; // re: NodeWrap

		
		
		
		template<typename NodeBaseT>
		template<class PyClass>
		void UnaryOpVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("set_child", &NodeBaseT::SetChild )
			.def("first_child", &NodeBaseT::first_child )
			;
		}
		
		template<typename NodeBaseT>
		template<class PyClass>
		void NaryOpVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("add_child", &NodeBaseT::AddChild )
			.def("first_child", &NodeBaseT::first_child )
			.def("children_size", &NodeBaseT::children_size )
			;
		}

		template<typename NodeBaseT>
		template<class PyClass>
		void SumMultOpVisitor<NodeBaseT>::visit(PyClass& cl) const
		{
			
			cl
			.def("add_child", addChild2)
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
			class_<Operator, boost::noncopyable, bases<Node>, std::shared_ptr<Operator> >("Operator", no_init)
			;

			// UnaryOperator class
			class_<UnaryOpWrap, boost::noncopyable, bases<Operator>, std::shared_ptr<UnaryOperator> >("UnaryOperator", no_init)
			.def(UnaryOpVisitor<UnaryOperator>())
			;

			// NaryOperator class
			class_<NaryOpWrap, boost::noncopyable, bases<Operator>, std::shared_ptr<NaryOperator> >("NaryOperator", no_init)
			.def(NaryOpVisitor<NaryOperator>())
			;

			// SumOperator class
			class_<SumOperator, bases<NaryOperator>, std::shared_ptr<SumOperator> >("SumOperator", init<const Nodeptr&, const Nodeptr &>() )
			.def(init<const Nodeptr&, bool, const Nodeptr&, bool>() )
			
			.def(SumMultOpVisitor<SumOperator>())
			;

			// NegateOperator class
			class_<NegateOperator, bases<UnaryOperator>, std::shared_ptr<NegateOperator> >("NegateOperator", init<const Nodeptr&>() )
			;

			// MultOperator class
			class_<MultOperator, bases<NaryOperator>, std::shared_ptr<MultOperator> >("MultOperator", init<const Nodeptr&, const Nodeptr &>() )
			.def(init<const Nodeptr&, bool, const Nodeptr&, bool>() )
			
			.def(SumMultOpVisitor<MultOperator>())
			;
			
			// PowerOperator class
			class_<PowerOperator, bases<Operator>, std::shared_ptr<PowerOperator> >("PowerOperator", init<const Nodeptr&, const Nodeptr &>() )
			
			.def(PowerOpVisitor<PowerOperator>())
			;
			
			// IntegerPowerOperator class
			class_<IntegerPowerOperator, bases<UnaryOperator>, std::shared_ptr<IntegerPowerOperator> >("IntegerPowerOperator", init<const Nodeptr&, int>()  )
			
			.def(PowerOpVisitor<IntegerPowerOperator>())
			;

			// SqrtOperator class
			class_<SqrtOperator, bases<UnaryOperator>, std::shared_ptr<SqrtOperator> >("SqrtOperator", init<const Nodeptr&>() )
			;

			// ExpOperator class
			class_<ExpOperator, bases<UnaryOperator>, std::shared_ptr<ExpOperator> >("ExpOperator", init<const Nodeptr&>() )
			;

			// LogOperator class
			class_<LogOperator, bases<UnaryOperator>, std::shared_ptr<LogOperator> >("LogOperator", init<const Nodeptr&>() )
			;
			
			

			// TrigOperator class
			class_<TrigOperator, boost::noncopyable, bases<Node>, std::shared_ptr<TrigOperator> >("TrigOperator", no_init)
			;

			// SinOperator class
			class_<SinOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<SinOperator> >("SinOperator", init<const Nodeptr&>() )
			;
			// CosOperator class
			class_<CosOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<CosOperator> >("CosOperator", init<const Nodeptr&>() )
			;
			// TanOperator class
			class_<TanOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<TanOperator> >("TanOperator", init<const Nodeptr&>() )
			;

			// ArcSinOperator class
			class_<ArcSinOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<ArcSinOperator> >("ArcSinOperator", init<const Nodeptr&>() )
			;
			// ArcCosOperator class
			class_<ArcCosOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<ArcCosOperator> >("ArcCosOperator", init<const Nodeptr&>() )
			;
			// ArcTanOperator class
			class_<ArcTanOperator, bases<TrigOperator, UnaryOperator>, std::shared_ptr<ArcTanOperator> >("ArcTanOperator", init<const Nodeptr&>() )
			;


			
			
		}
	} //namespace python
} // namespace bertini



