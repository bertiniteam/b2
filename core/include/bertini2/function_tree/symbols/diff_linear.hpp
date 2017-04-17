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
			 \brief Create a linear differential node with rational coefficients.  Only a linear, not a product.  Not homogenized, so no hom variable.
			 
			 \param coeffs_real Rational coefficients for the linear.
			 \param coeffs_imag Rational coefficients for the linear.
			 \param vars Variables associated with the linear.
			 
			 */
			DiffLinear(Mat<mpq_rational> const& coeffs_real, Mat<mpq_rational> const& coeffs_imag, VariableGroup const& vars) :
			coeffs_rat_real_(coeffs_real), coeffs_rat_imag_(coeffs_imag), variables_(vars), num_variables_(vars.size())
			{
				hom_variable_ = MakeInteger(0);
				is_homogenized_ = false;
			}

			
			/**
			 \brief Create a linear differential node with rational coefficients.  Only a linear, not a product.  Not homogenized, so no hom variable.
			 
			 \param coeffs_real Rational coefficients for the linear.
			 \param coeffs_imag Rational coefficients for the linear.
			 \param vars Variables associated with the linear.
			 
			*/
			DiffLinear(Mat<mpq_rational> const& coeffs_real, Mat<mpq_rational> const& coeffs_imag, VariableGroup const& vars) :
				coeffs_rat_real_(coeffs_real), coeffs_rat_imag_(coeffs_imag), variables_(vars), num_variables_(vars.size())
			{
				hom_variable_ = MakeInteger(0);
				is_homogenized_ = false;
			}

			/**
			 \brief Create a linear differential node with rational coefficients.  Only a linear, not a product.  With homogenizing variable.
			 
			 \param coeffs_real Rational coefficients for the linear.
			 \param coeffs_imag Rational coefficients for the linear.
			 \param vars Variables associated with the linear.
			 \param hom_var Homogenizing variable.
			 
			 */
			DiffLinear(Mat<mpq_rational> const& coeffs_real, Mat<mpq_rational> const& coeffs_imag, VariableGroup const& vars, std::shared_ptr<Variable> hom_var) :
			coeffs_rat_real_(coeffs_real), coeffs_rat_imag_(coeffs_imag), variables_(vars), hom_variable_(hom_var), num_variables_(vars.size())
			{
				is_homogenized_ = true;
			}

			
			/**
			 \brief Create a linear differential node with mpfr coefficients.  Only a linear, not a product.  Not homogenized, so no hom variable.
			 
			 \param coeffs mpfr coefficients for the linear.
			 \param vars Variables associated with the linear.
			 
			 */
			DiffLinear(Mat<mpfr> const& coeffs, VariableGroup const& vars) :
			variables_(vars), num_variables_(vars.size())
			{
				hom_variable_ = MakeInteger(0);
				is_homogenized_ = false;
				
				// Set the coefficient matrices with input matrices.
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				Mat<mpfr>& coeffs_mpfr_ref = std::get<Mat<mpfr>>(coeffs_);
				
				for(int ii = 0; ii < num_factors_; ++ii)
				{
					for(int jj = 0; jj < coeffs_mpfr.cols(); ++jj)
					{
						coeffs_dbl_ref(ii,jj) = static_cast<dbl>(coeffs_mpfr(ii,jj));
					}
				}
				coeffs_mpfr_ref = coeffs_mpfr;
			}
			
			/**
			 \brief Create a linear differential node with mpfr coefficients.  Only a linear, not a product.  With homogenizing variable.
			 
			 \param coeffs Rational coefficients for the linear.
			 \param vars Variables associated with the linear.
			 \param hom_var Homogenizing variable.
			 
			 */
			DiffLinear(Mat<mpfr> const& coeffs, VariableGroup const& vars, std::shared_ptr<Variable> hom_var) :
			variables_(vars), hom_variable_(hom_var), num_variables_(vars.size())
			{
				is_homogenized_ = true;
				
				// Set the coefficient matrices with input matrices.
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				Mat<mpfr>& coeffs_mpfr_ref = std::get<Mat<mpfr>>(coeffs_);

				for(int ii = 0; ii < num_factors_; ++ii)
				{
					for(int jj = 0; jj < coeffs_mpfr.cols(); ++jj)
					{
						coeffs_dbl_ref(ii,jj) = static_cast<dbl>(coeffs(ii,jj));
					}
				}
				coeffs_mpfr_ref = coeffs;

			}

			
			
			
			
			
			
			
			
			void SetupCoeffs(Mat<mpfr> const& coeffs)
			{
				// Set the coefficient matrices with input matrices.
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				Mat<mpfr>& coeffs_mpfr_ref = std::get<Mat<mpfr>>(coeffs_);
				
				for(int ii = 0; ii < num_factors_; ++ii)
				{
					for(int jj = 0; jj < coeffs_mpfr.cols(); ++jj)
					{
						coeffs_dbl_ref(ii,jj) = static_cast<dbl>(coeffs(ii,jj));
					}
				}
				coeffs_mpfr_ref = coeffs;
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
					for (int jj = 0; jj < num_variables_ + 1; ++jj)
					{
						std::cout << std::get< Mat<mpfr> >(coeffs_)(ii,jj) << " | ";
					}
					std::cout << "\n";
				}
			}
			
			
			
			
			
			
			
			
			
			/**
			 \brief Accept a variable and add it to the list of variables in linear factors.
			 
			 \param child Variable to be added.  If node is not a Variable, then runtime error is thrown.
			 
			 */
			//			void AddVariable(std::shared_ptr<Node> child) override
			//			{
			//				if(std::dynamic_pointer_cast< std::shared_ptr<Variable> >(child) == nullptr)
			//				{
			//					throw std::runtime_error("Only Variable node types can be children of a LinearProduct node.");
			//				}
			//			}
			
			
			
			
			
			
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
				std::shared_ptr<Node> n = MakeVariable("x");
				return n;
			}
			
			
			/**
			 \brief Computes the degree for a particular variable.  If that variable is part of the linear product, the degree is equal to the number of factors.
			 
			 \param v The variable we are determining the degree with respect to.
			 \return Degree of polynomial with respect to variable v.
			 */
			int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
			{
				int deg = 0;
				
				// If v is part of the linear product
				if(std::find(variables_.begin(), variables_.end(), v) != std::end(variables_))
				{
					deg = num_factors_;
				}
				
				return deg;
			};
			
			
			
			
			/**
			 \brief Computer the degree for a particular group of variables.  If one of the variables is a part of the linear product, the degree is equal to the number of factors.
			 
			 \param vars The group of variables we are determing the degree with respect to.
			 
			 \return Degree of polynomial with respect to variable group.
			 */
			
			int Degree(VariableGroup const& vars) const override
			{
				int deg = 0;
				
				for (auto v = vars.begin(); v != vars.end(); ++v)
				{
					// if v is a part of the linear product
					if(std::find(variables_.begin(), variables_.end(), *v) != std::end(variables_))
					{
						deg = num_factors_;
						break;
					}
				}
				
				return deg;
			};
			
			/**
			 \brief Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.
			 
			 \param vars The variable group computing degree with respect to.
			 \return Multidegree vector of polynomial with repect to variable group.
			 */
			std::vector<int> MultiDegree(VariableGroup const& vars) const override
			{
				std::vector<int> degs(vars.size(), 0);
				
				for (auto v = vars.begin(); v != vars.end(); ++v)
				{
					// If v is part of the linear product
					if(std::find(variables_.begin(), variables_.end(), *v) != std::end(variables_))
					{
						*(degs.begin()+(v-vars.begin())) = num_factors_;
					}
				}
				
				
				return degs;
			}
			
			
			
			
			
			/**
			 \brief Homogenize a sum, with respect to a variable group, and using a homogenizing variable.
			 
			 \param vars Variable group to homogenize with respect to.
			 \param homvar Homogenization variable.
			 */
			void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
			{
				//				std::cout << "Homogenizing linprod with num factors = " << num_factors_ << std::endl;  //DEBUGING
				
				// Check if vars is the variable group for this linear product
				bool is_vargroup_same = true; // Check if vars and variables_ are the same
				bool is_v_in_vars = false; // Check if at least one v in variables_ is in vars
				
				
				// Is homvar in vars?
				VariableGroup temp_vars = vars;
				temp_vars.erase(std::remove(temp_vars.begin(), temp_vars.end(), homvar), temp_vars.end());
				//				auto var_it = std::find(vars.begin(), vars.end(), homvar);
				//				// If yes, remove it from vars.
				//				if(var_it != std::end(vars))
				//				{
				//					std::cout << " it = " << **var_it << std::endl;
				////					temp_vars.erase(vars[0]);
				//					std::remove(temp_vars.begin(), temp_vars.end(), homvar);
				//				}
				//
				//				for (auto v : temp_vars)
				//				{
				//					std::cout << "v = " << *v << std::endl;
				//				}
				
				
				// Which group is larger, variables_ or vars?
				VariableGroup larger_group;
				VariableGroup smaller_group;
				if(variables_.size() > temp_vars.size())
				{
					larger_group = variables_;
					smaller_group = temp_vars;
				}
				else
				{
					larger_group = temp_vars;
					smaller_group = variables_;
				}
				
				// Is vars the same as variables_?
				for (auto const v : larger_group)
				{
					// If v is in vars?
					if(std::find(smaller_group.begin(), smaller_group.end(), v) != std::end(smaller_group))
					{
						is_v_in_vars = true;
					}
					else
					{
						is_vargroup_same = false;
					}
				}
				
				if(!is_vargroup_same && is_v_in_vars)
				{
					throw std::runtime_error("attempting to homogenize linear product with respect to part of the variables, but not all");
				}
				
				if(is_vargroup_same)
				{
					hom_variable_ = homvar;
					is_homogenized_ = true;
				}
				//				else
				//				{
				//					throw std::runtime_error("attemtping to homogenize linear product with respect to another variable group");
				//				}
				
			};
			
			bool IsHomogeneous(std::shared_ptr<Variable> const& v) const override
			{
				
				if(v == nullptr)
				{
					if(is_homogenized_)
					{
						return true;
					}
					else
					{
						return false;
					}
				}
				else
				{
					// If v is not part of the linear product
					if(std::find(variables_.begin(), variables_.end(), v) == std::end(variables_))
					{
						return true;
					}
					else
					{
						return false;
					}
				}
				
				
			};
			
			/**
				\brief Check for homogeneity, with respect to a variable group.  If vars is the same variable group as the one used to create the linear product, and it has been homogenized, then return true.
			 
				\param vars Variable group to check if homogenous with respect to.
			 
				\return boolean.
			 */
			bool IsHomogeneous(VariableGroup const& vars) const override
			{
				bool is_hom = false;
				
				// Check if vars is the variable group for this linear product
				bool is_vargroup_same = true; // Check if vars and variables_ are the same
				bool is_v_in_vars = false; // Check if at least one v in variables_ is in vars
				for (auto const v : variables_)
				{
					//					std::cout << *vars.begin() << std::endl;
					//					std::cout << *std::end( << std::endl;
					//					std::cout << *std::find(vars.begin(), vars.end(), variables_[0]) << std::endl;
					
					// If v is in vars?
					if(std::find(vars.begin(), vars.end(), v) != std::end(vars))
					{
						is_v_in_vars = true;
					}
					else
					{
						is_vargroup_same = false;
					}
				}
				
				if(is_homogenized_)
				{
					// If homvar is in vars?
					if(std::find(vars.begin(), vars.end(), hom_variable_) != std::end(vars))
					{
						is_v_in_vars = true;
					}
					else
					{
						is_vargroup_same = false;
					}
				}
				
				if(!is_vargroup_same && is_v_in_vars)
				{
					return false;
				}
				
				if(!is_vargroup_same)
				{
					is_hom = true;
				}
				else if(is_homogenized_)
				{
					is_hom = true;
				}
				
				return is_hom;
				
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
				evaluation_value = dbl(1);
				const Mat<dbl>& coeffs_ref = std::get< Mat<dbl> >(coeffs_);
				
				
				// Evaluate all afine variables and store
				for (int jj = 0; jj < num_variables_; ++jj)
				{
					variables_[jj]->EvalInPlace<dbl>(temp_var_d_[jj], diff_variable);
				}
				
				hom_variable_->EvalInPlace<dbl>(temp_var_d_[num_variables_], diff_variable);
				
				
				
				
				
				
				// Evaluate the linear product
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					// Add all terms in one linear factor and store in temp_sum_d_
					temp_sum_d_ = dbl(0);
					for (int jj = 0; jj < num_variables_ + 1; ++jj)
					{
						temp_sum_d_ += coeffs_ref(ii,jj)*temp_var_d_[jj];
					}// re: loop through variables
					
					
					// Multiply factors together
					evaluation_value *= temp_sum_d_;
				}// re: loop through factors
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
				evaluation_value.SetOne();
				const Mat<mpfr>& coeffs_ref = std::get< Mat<mpfr> >(coeffs_);
				
				
				// Evaluate all afine variables and store
				for (int jj = 0; jj < num_variables_; ++jj)
				{
					variables_[jj]->EvalInPlace<mpfr>(temp_var_mp_[jj], diff_variable);
				}
				hom_variable_->EvalInPlace<mpfr>(temp_var_mp_[num_variables_], diff_variable);
				
				
				
				
				
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					// Add all terms in one linear factor and store in temp_sum_d_
					temp_sum_mp_.SetZero();
					for (int jj = 0; jj < num_variables_ + 1; ++jj)
					{
						temp_sum_mp_ += coeffs_ref(ii,jj)*temp_var_mp_[jj];
					}
					
					
					// Multiply factors together
					evaluation_value *= temp_sum_mp_;
				}
			}
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		private:
			//			std::vector< std::vector< std::tuple<mpq_rational, dbl, mpfr> > > coeffs_;
			Mat<mpq_rational> coeffs_rat_real_;   ///< Matrix of real rational coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor, with the final column being the constant coefficient.  These rationals can then be downsampled for each data type.
			
			Mat<mpq_rational> coeffs_rat_imag_;   ///< Same as coeffs_rat_real_ but for imaginary portion of the coefficients.
			
			std::tuple< Mat<dbl>, Mat<mpfr> > coeffs_;   ///< Matrix of coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor with the final column being the constant coefficient.  This is a tuple with one matrix for each data type.
			
			VariableGroup variables_; ///< Variables to be used in each linear factor.  Does not have to correspond directly to a variable group from the system.
			
			std::shared_ptr<Node> hom_variable_; ///< The homogenizing variable for this variable group.  Initially this is set to an Integer = 1.  When we homogenize, this is set to the variable.
			
			//			std::vector< std::tuple< std::shared_ptr<Variable>, Vec<bool> > > hom_variables_; ///< A vector of tuples holding 1) the homogenizing variable and 2) which terms in the linear factor get multiplied by the hom variable.
			
			size_t num_variables_;  ///< The number of variables in each linear.
			bool is_rational_coeffs_; ///< Do we have a rational coefficient to downsample from?
			bool is_homogenized_;  ///< Have we homogenized the linear product?
			
			
			
			
			
			
			
			
			
		private:
			
			DiffLinear() = default;
			
			
			
			friend class boost::serialization::access;
			
			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<NaryOperator>(*this);
			}
			
			
			void PrecisionChangeSpecific(unsigned prec) const
			{
				temp_mp_.precision(prec);
			}
			
		};
		
		
		
	} // re: namespace node
} // re: namespace bertini




#endif
