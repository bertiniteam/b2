//This file is part of Bertini 2.
//
//python/system_export.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//python/system_export.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with python/system_export.hpp.  If not, see <http://www.gnu.org/licenses/>.
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
//  python/system_export.hpp:  Header file for exposing systems to python, including start systems.




#ifndef BERTINI_PYTHON_SYSTEM_EXPORT_HPP
#define BERTINI_PYTHON_SYSTEM_EXPORT_HPP
#include <bertini2/system/system.hpp>
#include <bertini2/system/start_systems.hpp>





#include "python_common.hpp"


namespace bertini{
	namespace python{
		
		using namespace bertini;
		
		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
		template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
		using VariableGroup = std::deque< std::shared_ptr<node::Variable> >;
		
		using dbl = std::complex<double>;
		using mpfr = bertini::complex;
		
		
		
		
		void ExportSystem();
		
		
		
		
		
		
		/**
		 System class 
		 */
		template<typename SystemBaseT>
		class SystemVisitor: public def_visitor<SystemVisitor<SystemBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:

			// precision functions
			void (bertini::System::*set_prec_)(unsigned) const = &bertini::System::precision;
			unsigned (bertini::System::*get_prec_)(void) const = &bertini::System::precision;

			void (bertini::System::*sysAddFunc1)(std::shared_ptr<node::Function> const&) = &bertini::System::AddFunction;
			void (bertini::System::*sysAddFunc2)(std::shared_ptr<node::Node> const&) = &bertini::System::AddFunction;
			

//			Vec<dbl> (System::*sysEval1)(const Vec<dbl> &) = &System::template Eval<dbl>;
//			Vec<dbl> (System::*sysEval1)(const Vec<dbl> &) = &System::template Eval<dbl>;
			
			std::vector<int> (bertini::System::*sysDeg1)() const = &bertini::System::Degrees;
			std::vector<int> (bertini::System::*sysDeg2)(VariableGroup const&) const = &bertini::System::Degrees;
			
			// Eval functions
			template <typename T>
			using Eval0_ptr = Vec<T> (SystemBaseT::*)() const;
			template <typename T>
			static Eval0_ptr<T> return_Eval0_ptr()
			{
				return &SystemBaseT::template Eval<T>;
			};

			template <typename T>
			using Eval1_ptr = Vec<T> (SystemBaseT::*)(const Vec<T>&) const;
			template <typename T>
			static Eval1_ptr<T> return_Eval1_ptr()
			{
				return &SystemBaseT::template Eval<T>;
			};
			
			template <typename T>
			using Eval2_ptr = Vec<T> (SystemBaseT::*)(const Vec<T>&, const T &) const;
			template <typename T>
			static Eval2_ptr<T> return_Eval2_ptr()
			{
				return &SystemBaseT::template Eval<T>;
			};
			
			
			// Jacobian Eval functions
			template <typename T>
			using Jac0_ptr = Mat<T> (SystemBaseT::*)() const;
			template <typename T>
			static Jac0_ptr<T> return_Jac0_ptr()
			{
				return &SystemBaseT::template Jacobian<T>;
			};

			template <typename T>
			using Jac1_ptr = Mat<T> (SystemBaseT::*)(const Vec<T>&) const;
			template <typename T>
			static Jac1_ptr<T> return_Jac1_ptr()
			{
				return &SystemBaseT::template Jacobian<T>;
			};
			
			template <typename T>
			using Jac2_ptr = Mat<T> (SystemBaseT::*)(const Vec<T>&, const T &) const;
			template <typename T>
			static Jac2_ptr<T> return_Jac2_ptr()
			{
				return &SystemBaseT::template Jacobian<T>;
			};

		};
		
		
		
		
		
		
		
		/**
		 StartSystem class 
		 */
		template<typename SystemBaseT>
		class StartSystemVisitor: public def_visitor<StartSystemVisitor<SystemBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
			
		private:
			template <typename T>
			using GenStart_ptr = Vec<T> (SystemBaseT::*)(unsigned long long) const;
			template <typename T>
			static GenStart_ptr<T> return_GenStart_ptr()
			{
				return &SystemBaseT::template StartPoint<T>;
			};
		};
		
	}
}


#endif
