//
//  system_export.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/9/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_system_export_hpp
#define Xcode_b2_system_export_hpp
#include <bertini2/system.hpp>





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
		
		
		
		
		///////// System class ////////////////
		template<typename SystemBaseT>
		class SystemVisitor: public def_visitor<SystemVisitor<SystemBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:
			void (bertini::System::*sysAddFunc1)(std::shared_ptr<node::Function> const&) = &bertini::System::AddFunction;
			void (bertini::System::*sysAddFunc2)(std::shared_ptr<node::Node> const&) = &bertini::System::AddFunction;
			

//			Vec<dbl> (System::*sysEval1)(const Vec<dbl> &) = &System::template Eval<dbl>;
//			Vec<dbl> (System::*sysEval1)(const Vec<dbl> &) = &System::template Eval<dbl>;
			
			std::vector<int> (bertini::System::*sysDeg1)() const = &bertini::System::Degrees;
			std::vector<int> (bertini::System::*sysDeg2)(VariableGroup const&) const = &bertini::System::Degrees;
			
			// Eval functions
			template <typename T>
			using Eval1_ptr = Vec<T> (SystemBaseT::*)(const Vec<T>&);
			template <typename T>
			static Eval1_ptr<T> return_Eval1_ptr()
			{
				return &SystemBaseT::template Eval<T>;
			};
			
			template <typename T>
			using Eval2_ptr = Vec<T> (SystemBaseT::*)(const Vec<T>&, const T &);
			template <typename T>
			static Eval2_ptr<T> return_Eval2_ptr()
			{
				return &SystemBaseT::template Eval<T>;
			};
			
			
			// Jacobian Eval functions
			template <typename T>
			using Jac1_ptr = Mat<T> (SystemBaseT::*)(const Vec<T>&);
			template <typename T>
			static Jac1_ptr<T> return_Jac1_ptr()
			{
				return &SystemBaseT::template Jacobian<T>;
			};
			
			template <typename T>
			using Jac2_ptr = Mat<T> (SystemBaseT::*)(const Vec<T>&, const T &);
			template <typename T>
			static Jac2_ptr<T> return_Jac2_ptr()
			{
				return &SystemBaseT::template Jacobian<T>;
			};

			


		};
		
	}
}


#endif
