//This file is part of Bertini 2.
//
//variable.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//diff_linear.hpp is distributed in the hope that it will be useful,
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
//  Created by Collins, James B. on 4/17/2017.
//
//
// diff_linear.hpp:  Declares the class DiffLinear.


/**
 \file linear_product.hpp
 
 \brief Provides the LinearProduct Node class.
 
 */
#ifndef BERTINI_FUNCTION_TREE_DIFFLIN_HPP
#define BERTINI_FUNCTION_TREE_DIFFLIN_HPP

#include "bertini2/function_tree/symbols/symbol.hpp"
#include "bertini2/function_tree/symbols/linear_product.hpp"
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
		class DiffLinear : public virtual Symbol
		{
		public:
			virtual ~DiffLinear() = default;
			
			
			/**
			 \brief Create a linear differential node. Only a linear, not a product.
			 
			 \param linear The linear that we are differentiating.
			 
			 */
			DiffLinear(std::shared_ptr<LinearProduct> const& linear)
			{
				// Set differentials of all variables in the linear
				linear->GetVariables(variables_);

				linear->GetHomVariable(hom_variable_);
				
				
				// Set coefficients with rational or mpfr type
				size_t num_variables_ = variables_.size();
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				Mat<mpfr>& coeffs_mpfr_ref = std::get<Mat<mpfr>>(coeffs_);
				coeffs_dbl_ref.resize(1, num_variables_+1);
				coeffs_mpfr_ref.resize(1, num_variables_+1);

				if(linear->IsRationalCoefficients())
				{
					coeffs_rat_real_.resize(1, num_variables_+1);
					coeffs_rat_real_.resize(1, num_variables_+1);
					linear->GetRatCoeffs(coeffs_rat_real_, coeffs_rat_imag_);
					for(int jj = 0; jj < num_variables_+1; ++jj)
					{
						coeffs_dbl_ref(0,jj).real( static_cast<double>(coeffs_rat_real_(0,jj)) );
						coeffs_dbl_ref(0,jj).imag( static_cast<double>(coeffs_rat_imag_(0,jj)) );
						coeffs_mpfr_ref(0,jj).real( static_cast<mpfr_float>(coeffs_rat_real_(0,jj)) );
						coeffs_mpfr_ref(0,jj).imag( static_cast<mpfr_float>(coeffs_rat_imag_(0,jj)) );
					}
				}
				else
				{
					linear->GetMPFRCoeffs(coeffs_mpfr_ref);
					
					for(int jj = 0; jj < num_variables_+1; ++jj)
					{
						coeffs_dbl_ref(jj) = static_cast<dbl>(coeffs_mpfr_ref(jj));
					}

				}
				
				
			}

			
			
			
			
			
			
			
			
			
			
			
			
			
			/**
			 \brief Reset variable values in this node
			 */
			void Reset() const override {};
			
			
			/**
			 Method for printing to output stream
			 */
			void print(std::ostream & target) const override{};
			
			
			
			
			/**
			 Return SumOperator whose children are derivatives of children_
			 */
			std::shared_ptr<Node> Differentiate() const override
			{
				return MakeInteger(0);
			}
			
			
			/**
			 \brief Computes the degree for a particular variable.  If that variable is part of the linear product, the degree is equal to the number of factors.
			 
			 \param v The variable we are determining the degree with respect to.
			 \return Degree of polynomial with respect to variable v.
			 */
			int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
			{
				return 0;
			};
			
			
			
			
			/**
			 \brief Computer the degree for a particular group of variables.  If one of the variables is a part of the linear product, the degree is equal to the number of factors.
			 
			 \param vars The group of variables we are determing the degree with respect to.
			 
			 \return Degree of polynomial with respect to variable group.
			 */
			
			int Degree(VariableGroup const& vars) const override
			{
				return 0;
			};
			
			/**
			 \brief Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.
			 
			 \param vars The variable group computing degree with respect to.
			 \return Multidegree vector of polynomial with repect to variable group.
			 */
			std::vector<int> MultiDegree(VariableGroup const& vars) const override
			{
				return std::vector<int>(vars.size(),0);
			}
			
			
			
			
			
			/**
			 \brief Homogenize a sum, with respect to a variable group, and using a homogenizing variable.
			 
			 \param vars Variable group to homogenize with respect to.
			 \param homvar Homogenization variable.
			 */
			void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
			{
			};
			
			/**
				\brief Check for homogeneity, with respect to a variable group.  If vars is the same variable group as the one used to create the linear product, and it has been homogenized, then return true.
			 
				\param vars Variable group to check if homogenous with respect to.
			 
				\return boolean.
			 */
			bool IsHomogeneous(VariableGroup const& vars) const override
			{
				return true;
			};
			
			
			
			/**
			 Change the precision of this variable-precision tree node.
			 
			 \param prec the number of digits to change precision to.
			 */
			virtual void precision(unsigned int prec) const {};
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		protected:
			/**
			 \brief Evaluation of linear product node.  Returns evaluation value.
			 
			 \param diff_variable Variable that we are differentiating with respect to.  Only for evaluating Jacobians.
			 */
			dbl FreshEval_d(std::shared_ptr<Variable> const& diff_variable) const override
			{
				dbl eval_value;
				
				this->FreshEval_d(eval_value, diff_variable);
				return eval_value;
			}
			
			/**
			 \brief Evaluation of linear product node IN PLACE.  Returns evaluation value.
			 
			 \param evaluation_value The in place variable that stores the evaluation.
			 \param diff_variable Variable that we are differentiating with respect to.  Only for evaluating Jacobians.
			 */
			void FreshEval_d(dbl& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override
			{
				for(int ii = 0; ii < variables_.size(); ++ii)
				{
					if(diff_variable == variables_[ii])
					{
						auto& coeff_ref = std::get<Mat<dbl>>(coeffs_);
						evaluation_value = coeff_ref(0,ii);
						return;
					}
				}
				
				// If not one of the affine variables
				if(diff_variable == hom_variable_)
				{
					auto& coeff_ref = std::get<Mat<dbl>>(coeffs_);
					evaluation_value = coeff_ref(0,variables_.size()-1);
					return;
				}
				
				
				// If none of the variables
				evaluation_value = dbl(0);
				return;
			}
			
			
			/**
			 \brief Evaluation of linear product node.  Returns evaluation value.
			 
			 \param diff_variable Variable that we are differentiating with respect to.  Only for evaluating Jacobians.
			 */
			mpfr FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override
			{
				mpfr eval_value;
				
				this->FreshEval_mp(eval_value, diff_variable);
				return eval_value;
			}
			
			
			/**
			 \brief Evaluation of linear product node IN PLACE.  Returns evaluation value.
			 
			 \param evaluation_value The in place variable that stores the evaluation.
			 \param diff_variable Variable that we are differentiating with respect to.  Only for evaluating Jacobians.
			 */
			void FreshEval_mp(mpfr& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override
			{
				for(int ii = 0; ii < variables_.size(); ++ii)
				{
					if(diff_variable == variables_[ii])
					{
						auto& coeff_ref = std::get<Mat<mpfr>>(coeffs_);
						evaluation_value = coeff_ref(0,ii);
						return;
					}
				}
				
				// If not one of the affine variables
				if(diff_variable == hom_variable_)
				{
					auto& coeff_ref = std::get<Mat<mpfr>>(coeffs_);
					evaluation_value = coeff_ref(0,variables_.size()-1);
					return;
				}
				
				
				// If none of the variables
				evaluation_value.SetZero();
				return;
			}
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		private:
			//			std::vector< std::vector< std::tuple<mpq_rational, dbl, mpfr> > > coeffs_;
			Mat<mpq_rational> coeffs_rat_real_;   ///< Matrix of real rational coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor, with the final column being the constant coefficient.  These rationals can then be downsampled for each data type.
			
			Mat<mpq_rational> coeffs_rat_imag_;   ///< Same as coeffs_rat_real_ but for imaginary portion of the coefficients.
			
			std::tuple< Mat<dbl>, Mat<mpfr> > coeffs_;   ///< Matrix of coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor with the final column being the constant coefficient.  This is a tuple with one matrix for each data type.
			
			VariableGroup variables_; ///< Differentials of variables used in the linear.
			
			std::shared_ptr<Node> hom_variable_; ///< The homogenizing variable for this variable group.  Initially this is set to an Integer = 0.  When we homogenize, this is set to the variable.
			
			
			size_t num_variables_;  ///< The number of variables in each linear.
			bool is_rational_coeffs_; ///< Do we have a rational coefficient to downsample from?
			
			
			
			
			
			
			
			
			
		private:
			
			DiffLinear() = default;
			
			
			
			friend class boost::serialization::access;
			
			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<NaryOperator>(*this);
			}
			
			
			void PrecisionChangeSpecific(unsigned prec) const
			{
			}
			
		};
		
		
		
	} // re: namespace node
} // re: namespace bertini




#endif
