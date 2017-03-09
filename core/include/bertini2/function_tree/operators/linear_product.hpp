//This file is part of Bertini 2.
//
//variable.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//variable.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with variable.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2015
//
//
//  Created by Collins, James B. on 3/6/2017.
//
//
// linear_product.hpp:  Declares the class LinearProduct.


/**
 \file linear_product.hpp
 
 \brief Provides the LinearProduct Node class.
 
 */
#ifndef BERTINI_FUNCTION_TREE_LINPRODUCT_HPP
#define BERTINI_FUNCTION_TREE_LINPRODUCT_HPP

#include "bertini2/function_tree/symbols/symbol.hpp"
#include "bertini2/function_tree/factory.hpp"
#include "bertini2/eigen_extensions.hpp"

template<typename NumType> using Mat = bertini::Mat<NumType>;

namespace  bertini {
	namespace node{
		/**
		 \brief Represents a product of linears as leaves in the function tree.
		 
		 This class represents a product of linear factors of a fixed set of variables.  This could
		 also be a single linear of a fixed set of variables.
		 
		 When differentiated, produces a differential referring to it.
		 */
		class LinearProduct : public virtual NaryOperator
		{
		public:
			virtual ~LinearProduct() = default;
			
			
			
			/**
			 /brief Constructor for a linear product node that generates random coefficients automatically.
			 
			 \param variables A deque of variable nodes that are used in each linear factor.  This does not have to be an actual system variable group.
			 \param num_factors The number of linear factors in the product.
			 
			*/
			LinearProduct(VariableGroup variables, int num_factors) : variables_(variables), num_factors_(num_factors)
			{
				num_variables_ = variables.size();
				
				// Resize coeffs matrix
				coeffs_rat_real_.resize(num_factors_, num_variables_);
				coeffs_rat_imag_.resize(num_factors_, num_variables_);
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				coeffs_dbl_ref.resize(num_factors_, num_variables_);
				Mat<mpfr>& coeffs_mpfr_ref = std::get<Mat<mpfr>>(coeffs_);
				coeffs_mpfr_ref.resize(num_factors_, num_variables_);
				
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					for (int jj = 0; jj < num_variables_; ++jj)
					{
						// Generate random constants as mpq_rationals.  Then downsample to mpfr and dbl.
						// TODO: RandomRat() does not generate random numbers.  Same each run.
						coeffs_rat_real_(ii,jj) = RandomRat();
						coeffs_rat_imag_(ii,jj) = RandomRat();
						coeffs_dbl_ref(ii,jj).real( static_cast<double>(coeffs_rat_real_(ii,jj)) );
						coeffs_dbl_ref(ii,jj).imag( static_cast<double>(coeffs_rat_imag_(ii,jj)) );
						coeffs_mpfr_ref(ii,jj).real( static_cast<mpfr_float>(coeffs_rat_real_(ii,jj)) );
						coeffs_mpfr_ref(ii,jj).imag( static_cast<mpfr_float>(coeffs_rat_imag_(ii,jj)) );
					}
					
				}
			}
			

			
			/**
			 /brief Constructor for a linear product node that generates random coefficients automatically.
			 
			 \param variables A deque of variable nodes that are used in each linear factor.  This does not have to be an actual system variable group.
			 \param num_factors The number of linear factors in the product.
			 \param coeffs_real A matrix of dbl data types for all the coefficients in the linear factors.
			 \param coeffs_mpfr A matrix of mpfr data types for all the coefficients in the linear factors.
			 
			 */
			LinearProduct(VariableGroup variables, int num_factors, Mat<dbl>& coeffs_dbl, Mat<mpfr>& coeffs_mpfr)
					: variables_(variables), num_factors_(num_factors)
			{
				num_variables_ = variables.size();
				
				// Resize coeffs matrix
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				coeffs_dbl_ref.resize(num_factors_, num_variables_);
				Mat<mpfr>& coeffs_mpfr_ref = std::get<Mat<mpfr>>(coeffs_);
				coeffs_mpfr_ref.resize(num_factors_, num_variables_);
				
				
				// Set the coefficient matrices with input matrices.
				coeffs_dbl_ref = coeffs_dbl;
				coeffs_mpfr_ref = coeffs_mpfr;
				
			}

			
			//////////////////////////////////////////
			//
			//         Testing/Debugging
			//
			//////////////////////////////////////////
			void print_coeffs()
			{
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					for (int jj = 0; jj < num_variables_; ++jj)
					{
						std::cout << std::get< Mat<mpfr> >(coeffs_)(ii,jj) << " | ";
					}
					std::cout << "\n";
				}
			}
			
			
			
			
//			SumOperator(const std::shared_ptr<Node> & s, bool add_or_sub)
//			{
//				AddChild(s, add_or_sub);
//			}
//			
//			SumOperator(const std::shared_ptr<Node> & left, const std::shared_ptr<Node> & right)
//			{
//				AddChild(left);
//				AddChild(right);
//			}
//			
//			
//			SumOperator(const std::shared_ptr<Node> & left, bool add_or_sub_left, const std::shared_ptr<Node> & right, bool add_or_sub_right)
//			{
//				AddChild(left, add_or_sub_left);
//				AddChild(right, add_or_sub_right);
//			}
//			
//			
//			SumOperator& operator+=(const std::shared_ptr<Node> & rhs)
//			{
//				this->AddChild(rhs);
//				return *this;
//			}
//			
//			SumOperator& operator-=(const std::shared_ptr<Node> & rhs)
//			{
//				this->AddChild(rhs,false);
//				return *this;
//			}
			
			
			
			
			
			
			//Special Behaviour: by default all terms added are positive
			void AddChild(std::shared_ptr<Node> child) override
			{
				
			}
			
			
			
			//Special Behaviour: Pass bool to set sign of term: true = add, false = subtract
			void AddChild(std::shared_ptr<Node> child, bool sign) // not an override
			{
				
			}
			
			
			/**
			 Method for printing to output stream
			 */
			void print(std::ostream & target) const override{};
			
			
			
			
			/**
			 Return SumOperator whose children are derivatives of children_
			 */
			std::shared_ptr<Node> Differentiate() const override
			{
				std::shared_ptr<Node> n = MakeVariable("x");
				return n;
			}
			
			
			/**
			 Compute the degree of a node.  For sum functions, the degree is the max among summands.
			 */
			int Degree(std::shared_ptr<Variable> const& v = nullptr) const override{return 0;};
			
			
			int Degree(VariableGroup const& vars) const override{return 0;};
			
			/**
			 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.
			 */
			std::vector<int> MultiDegree(VariableGroup const& vars) const override
			{
				std::vector<int> v(1);
				return v;
			}
			
			
			
			
			
			/**
			 Homogenize a sum, with respect to a variable group, and using a homogenizing variable.
			 */
			void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override{};
			
			bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override{return false;};
			
			/**
			 Check for homogeneity, with respect to a variable group.
			 */
			bool IsHomogeneous(VariableGroup const& vars) const override {return false;};
			
			
			
			
			
			
		protected:
			/**
			 Specific implementation of FreshEval for add and subtract.
			 If child_sign_ = true, then add, else subtract
			 */
			dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override{return dbl(2);};
			
			/**
			 Specific implementation of FreshEval in place for add and subtract.
			 If child_sign_ = true, then add, else subtract
			 */
			void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override{};
			
			
			/**
			 Specific implementation of FreshEval for add and subtract.
			 If child_sign_ = true, then add, else subtract
			 */
			mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override{return mpfr(2);};
			
			/**
			 Specific implementation of FreshEval for add and subtract.
			 If child_sign_ = true, then add, else subtract
			 */
			void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override{};
			
		private:
//			std::vector< std::vector< std::tuple<mpq_rational, dbl, mpfr> > > coeffs_;
			Mat<mpq_rational> coeffs_rat_real_;   ///< Matrix of rational coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor.  These rationals can then be downsampled for each data type.
			Mat<mpq_rational> coeffs_rat_imag_;   ///< Same as coeffs_rat_real_ but for imaginary portion of the coefficients.
			std::tuple< Mat<dbl>, Mat<mpfr> > coeffs_;   ///< Matrix of coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor.  This is a tuple with one matrix for each data type.
			VariableGroup variables_; ///< Variables to be used in each linear factor.  Does not have to correspond directly to a variable group from the system.
			
			size_t num_factors_;  ///< The number of factors in the linear product.
			size_t num_variables_;  ///< The number of variables in each linear.
			
			
			
			
			
			
			
			
			
		private:
			
			LinearProduct() = default;
			
			friend class boost::serialization::access;
			
			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<NaryOperator>(*this);
			}
			
			
			void PrecisionChangeSpecific(unsigned prec) const override
			{
				temp_mp_.precision(prec);
			}
			
			mutable mpfr temp_mp_;
			mutable dbl temp_d_;
		};
		
		
		
	} // re: namespace node
} // re: namespace bertini




#endif
