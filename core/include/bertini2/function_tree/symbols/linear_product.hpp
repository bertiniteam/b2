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
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
//  James Collins
//  West Texas A&M University
//  Spring, Summer 2017
//
//
// silviana amethyst
// University of Wisconsin - Eau Claire
//
//  Created by Collins, James B. on 3/6/2017.


/**
 \file linear_product.hpp
 
 \brief Provides the LinearProduct and DiffLinear Node classes.
 
 */
#ifndef BERTINI_FUNCTION_TREE_LINPRODUCT_HPP
#define BERTINI_FUNCTION_TREE_LINPRODUCT_HPP

#include "bertini2/function_tree.hpp"

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
		class LinearProduct : public virtual Symbol, public virtual EnableSharedFromThisVirtual<LinearProduct>
		{
		public:
			BERTINI_DEFAULT_VISITABLE()

			virtual ~LinearProduct() = default;
			
			
			template<typename... Ts> 
			static 
			std::shared_ptr<LinearProduct> Make(Ts&& ...ts){ 
				return std::shared_ptr<LinearProduct>( new LinearProduct(ts...) );
			}


		private:
			/**
			 /brief Constructor for a linear product node that generates random coefficients automatically.
			 
			 \param variables A deque of variable nodes that are used in each linear factor.  This does not have to be an actual system variable group.
			 \param num_factors The number of linear factors in the product.
			 
			 */
			LinearProduct(VariableGroup variables, int num_factors, bool is_hom_vars = false) :
				variables_(variables), num_factors_(num_factors), is_hom_vars_(is_hom_vars)
			{
				//				num_variables_ = variables.size();  //Include homogenize variable
				//
				//				// Resize coeffs matrix
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				//				coeffs_dbl_ref.resize(num_factors_, num_variables_+1);
				Mat<mpfr_complex>& coeffs_mpfr_ref = std::get<Mat<mpfr_complex>>(coeffs_);
				//				coeffs_mpfr_ref.resize(num_factors_, num_variables_+1);
				//				// Resize temporary variable holders and fill constant with 1
				//				temp_var_d_.resize(num_variables_ + 1);
				//				temp_var_d_[num_variables_] = dbl(1);
				//				temp_var_mp_.resize(num_variables_ + 1);
				//				temp_var_mp_[num_variables_] = mpfr_complex(1);
				
				SetupVariables(num_factors, variables);
				coeffs_rat_real_.resize(num_factors_, num_variables_+1);
				coeffs_rat_imag_.resize(num_factors_, num_variables_+1);
				
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					for (int jj = 0; jj < num_variables_+1; ++jj)
					{
						// Generate random constants as mpq_rationals.  Then downsample to mpfr_complex and dbl.
						// TODO: RandomRat() does not generate random numbers.  Same each run.
						coeffs_rat_real_(ii,jj) = RandomRat();
						coeffs_rat_imag_(ii,jj) = RandomRat();
						coeffs_dbl_ref(ii,jj).real( static_cast<double>(coeffs_rat_real_(ii,jj)) );
						coeffs_dbl_ref(ii,jj).imag( static_cast<double>(coeffs_rat_imag_(ii,jj)) );
						coeffs_mpfr_ref(ii,jj).real( static_cast<mpfr_float>(coeffs_rat_real_(ii,jj)) );
						coeffs_mpfr_ref(ii,jj).imag( static_cast<mpfr_float>(coeffs_rat_imag_(ii,jj)) );
					}
				}
				
				is_rational_coeffs_ = true;
				
                // If we are creating a linear product for variable group that is already homogenized,
                //   we do not use the constant(last) column in the coefficient matrix
				if(is_hom_vars)
				{
					is_homogenized_ = true;
					hom_variable_ = Integer::Make(0);
					for(int ii = 0; ii < coeffs_mpfr_ref.rows(); ++ii)
					{
						coeffs_rat_real_(ii,num_variables_) = mpq_rational(0);
						coeffs_rat_imag_(ii,num_variables_) = mpq_rational(0);
						coeffs_dbl_ref(ii,num_variables_) = dbl(0);
						coeffs_mpfr_ref(ii,num_variables_) = mpfr_complex(0);
					}
				}

				
			}
			
			
			
			/**
			 \brief Constructor for a linear product node that passes in random coefficients.
			 
			 \param variables A deque of variable nodes that are used in each linear factor.  This does not have to be an actual system variable group.
			 \param num_factors The number of linear factors in the product.
			 \param coeffs_real A matrix of dbl data types for all the coefficients in the linear factors.  Matrix must be of size num_factors x (num_variables + 1) with constant coefficient in last column.
			 \param coeffs_mpfr A matrix of mpfr_complex data types for all the coefficients in the linear factors.  Matrix must be of size num_factors x (num_variables + 1) with constant coefficient in last column.
			 
			 */
			LinearProduct(VariableGroup const& variables, Mat<mpfr_complex> const& coeffs_mpfr, bool is_hom_vars = false)
			: variables_(variables), num_factors_(coeffs_mpfr.rows()), is_hom_vars_(is_hom_vars)
			{
				if(variables.size()+1 != coeffs_mpfr.cols())
					throw std::runtime_error("attempting to construct a linear product manually, not enough columns for the number of variables in the variable group");
				
				
				// Resize coeffs matrix
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				Mat<mpfr_complex>& coeffs_mpfr_ref = std::get<Mat<mpfr_complex>>(coeffs_);
				
				SetupVariables(num_factors_, variables);
				
				coeffs_mpfr_ref = coeffs_mpfr;
                
                // If we are creating a linear product for variable group that is already homogenized,
                //   we do not use the constant(last) column in the coefficient matrix
				if(is_hom_vars)
				{
					is_homogenized_ = true;
					hom_variable_ = Integer::Make(0);
					for(int ii = 0; ii < coeffs_mpfr_ref.rows(); ++ii)
					{
						coeffs_mpfr_ref(ii,num_variables_) = mpfr_complex(0);
					}
				}
				
				// Set the coefficient matrices with input matrices.
				
				for(int ii = 0; ii < num_factors_; ++ii)
				{
					for(int jj = 0; jj < coeffs_mpfr.cols(); ++jj)
					{
						coeffs_dbl_ref(ii,jj) = static_cast<dbl>(coeffs_mpfr_ref(ii,jj));
					}
				}
				
				
				is_rational_coeffs_ = false;
				
			}
			
		

		public:
			
			
			
			
			
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
			void Reset() const override
			{
				Node::ResetStoredValues();
				for(auto& v : variables_)
				{
					v->Reset();
				}
			};
			
			
			/**
			 Method for printing to output stream
			 */
			void print(std::ostream & target) const override;
			
			
			
			
			/**
			 Return SumOperator whose children are derivatives of children_
			 */
			std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override    ;
	
			
			
			
			/**
			 \brief Computes the degree for a particular variable.  If that variable is part of the linear product, the degree is equal to the number of factors.
			 
			 \param v The variable we are determining the degree with respect to.
			 \return Degree of polynomial with respect to variable v.
			 */
			int Degree(std::shared_ptr<Variable> const& v = nullptr) const override;
			
			
			
			
			/**
			 \brief Computer the degree for a particular group of variables.  If one of the variables is a part of the linear product, the degree is equal to the number of factors.
			 
			 \param vars The group of variables we are determing the degree with respect to.
			 
			 \return Degree of polynomial with respect to variable group.
			 */
			
			int Degree(VariableGroup const& vars) const override;
			
			
			/**
			 \brief Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.
			 
			 \param vars The variable group computing degree with respect to.
			 \return Multidegree vector of polynomial with repect to variable group.
			 */
			std::vector<int> MultiDegree(VariableGroup const& vars) const override;
			
			
			
			
			
			
			
			/**
			 \brief Homogenize a sum, with respect to a variable group, and using a homogenizing variable.
			 
			 \param vars Variable group to homogenize with respect to.
			 \param homvar Homogenization variable.
			 */
			void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override;
			
			
			
			/**
				\brief Check for homogeneity, with respect to a variable.  If vars is the same variable group as the one used to create the linear product, and it has been homogenized, then return true.
			 
				\param vars Variable group to check if homogenous with respect to.
			 
				\return boolean.
			 */

			bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override;
			
			/**
				\brief Check for homogeneity, with respect to a variable group.  If vars is the same variable group as the one used to create the linear product, and it has been homogenized, then return true.
			 
				\param vars Variable group to check if homogenous with respect to.
			 
				\return boolean.
			 */
			bool IsHomogeneous(VariableGroup const& vars) const override;
			
			
			
			/**
			 Change the precision of this variable-precision tree node.
			 
			 \param prec the number of digits to change precision to.
			 */
			virtual void precision(unsigned int prec) const override
			{
				auto& val_pair = std::get< std::pair<mpfr_complex,bool> >(current_value_);
				val_pair.first.precision(prec);
				
				this->PrecisionChangeSpecific(prec);
				
				for (auto iter : variables_)
					iter->precision(prec);

			};
			
			

            unsigned EliminateZeros() override { return 0;};
            
            unsigned EliminateOnes() override {return 0;};
			
			
			
			
			
			
			
			
			
			
			
			
			

			
			
			
			/**
			 \brief Getter for the differential of all variables
			 
			*/
			void GetVariables(VariableGroup& vars) const
			{
				vars = variables_;
			}
			
			
			
			
			/**
			 \brief Getter for the homogenizing variable
			*/
			void GetHomVariable(std::shared_ptr<Node>& hom_var) const
			{
				hom_var = hom_variable_;
			}
			
			
			
			/**
			 \brief Getter for rational coefficients
			*/
			void GetRatCoeffs(Mat<mpq_rational>& coeffs_real, Mat<mpq_rational>& coeffs_imag) const
			{
				coeffs_real = coeffs_rat_real_;
				coeffs_imag = coeffs_rat_imag_;
			}
			
			/**
			 \brief Getter for mpfr_complex coefficients
			*/
			void GetMPFRCoeffs(Mat<mpfr_complex>& coeffs) const
			{
				auto& coeffs_mpfr_ref = std::get<Mat<mpfr_complex>>(coeffs_);
				coeffs = coeffs_mpfr_ref;
			}
			
			/**
			 \brief Are there rational coefficients?
			*/
			bool IsRationalCoefficients()
			{
				return is_rational_coeffs_;
			}
			
			
			
			/**
			 \brief Getter for the coefficients of a single linear factor.  Constant coefficient is last element of the Vector.
			 
			 \param index Index of the linear factor.
			 
			 \return Vector containing coefficients.
			 */
			template<typename CType>
			Vec<CType> GetCoeffs(size_t index)
			{
				Vec<CType> coeff_ret(num_variables_ + 1);
				Mat<CType>& coeff_ref = std::get<Mat<CType>>(coeffs_);
				
				for(int jj = 0; jj < num_variables_ + 1; ++jj)
				{
					coeff_ret(jj) = coeff_ref(index,jj);
				}
				
				return coeff_ret;
			}
			
			


			
			

			
			
			
			
			
			
			
			
			
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
			mpfr_complex FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override
			{
				mpfr_complex eval_value;
				
				this->FreshEval_mp(eval_value, diff_variable);
				return eval_value;
			}
			
			
			/**
			 \brief Evaluation of linear product node IN PLACE.  Returns evaluation value.
			 
			 \param evaluation_value The in place variable that stores the evaluation.
			 \param diff_variable Variable that we are differentiating with respect to.  Only for evaluating Jacobians.
			 */
			void FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override
			{
				evaluation_value = mpfr_complex(1);
				const Mat<mpfr_complex>& coeffs_ref = std::get< Mat<mpfr_complex> >(coeffs_);
				
				// Evaluate all afine variables and store
				for (int jj = 0; jj < num_variables_; ++jj)
				{
					variables_[jj]->EvalInPlace<mpfr_complex>(temp_var_mp_[jj], diff_variable);
				}
				hom_variable_->EvalInPlace<mpfr_complex>(temp_var_mp_[num_variables_], diff_variable);
				
				
				
				
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					// Add all terms in one linear factor and store in temp_sum_d_
					temp_sum_mp_ = mpfr_complex(0);
					for (int jj = 0; jj < num_variables_ + 1; ++jj)
					{
						temp_sum_mp_ += coeffs_ref(ii,jj)*temp_var_mp_[jj];
					}
					
					
					// Multiply factors together
					evaluation_value *= temp_sum_mp_;
				}
			}
			
			
		
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		private:
			Mat<mpq_rational> coeffs_rat_real_;   ///< Matrix of real rational coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor, with the final column being the constant coefficient.  These rationals can then be downsampled for each data type.
			
			Mat<mpq_rational> coeffs_rat_imag_;   ///< Same as coeffs_rat_real_ but for imaginary portion of the coefficients.
			
			std::tuple< Mat<dbl>, Mat<mpfr_complex> > coeffs_;   ///< Matrix of coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor with the final column being the constant coefficient.  This is a tuple with one matrix for each data type.
			
			VariableGroup variables_; ///< Variables to be used in each linear factor.  Does not have to correspond directly to a variable group from the system.
			
			std::shared_ptr<Node> hom_variable_; ///< The homogenizing variable for this variable group.  Initially this is set to an Integer = 1.  When we homogenize, this is set to the variable.
			
			
			size_t num_factors_;  ///< The number of factors in the linear product.
			size_t num_variables_;  ///< The number of variables in each linear.
			bool is_rational_coeffs_; ///< Do we have a rational coefficient to downsample from?
			bool is_homogenized_ = false;  ///< Have we homogenized the linear product?
			bool is_hom_vars_ = false;  ///< Is this linear for a homogeneous variable group?
			
			
			mutable std::vector<dbl> temp_var_d_;
			mutable std::vector<mpfr_complex> temp_var_mp_;
			mutable mpfr_complex temp_sum_mp_;
			mutable dbl temp_sum_d_;
			
			
			
			
			
			
			
		private:
			
			LinearProduct() = default;
			
			LinearProduct(VariableGroup const& variables, std::shared_ptr<Node> const& hom_var,
						  Mat<mpq_rational> const& coeffs_real, Mat<mpq_rational> const& coeffs_imag, bool is_hom_vars) :
				variables_(variables), num_factors_(coeffs_real.rows()), hom_variable_(hom_var), is_hom_vars_(is_hom_vars)
			{
				num_variables_ = variables.size();
				
				// Resize coeffs matrix
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				coeffs_dbl_ref.resize(num_factors_, num_variables_+1);
				Mat<mpfr_complex>& coeffs_mpfr_ref = std::get<Mat<mpfr_complex>>(coeffs_);
				coeffs_mpfr_ref.resize(num_factors_, num_variables_+1);
				
				// Resize temporary variable holders
				temp_var_d_.resize(num_variables_ + 1);
				temp_var_d_[num_variables_] = dbl(1);
				temp_var_mp_.resize(num_variables_ + 1);
				temp_var_mp_[num_variables_] = mpfr_complex(1);
				
				
				coeffs_rat_real_ = coeffs_real;
				coeffs_rat_imag_ = coeffs_imag;
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					for (int jj = 0; jj < num_variables_+1; ++jj)
					{
						coeffs_dbl_ref(ii,jj).real( static_cast<double>(coeffs_rat_real_(ii,jj)) );
						coeffs_dbl_ref(ii,jj).imag( static_cast<double>(coeffs_rat_imag_(ii,jj)) );
						coeffs_mpfr_ref(ii,jj).real( static_cast<mpfr_float>(coeffs_rat_real_(ii,jj)) );
						coeffs_mpfr_ref(ii,jj).imag( static_cast<mpfr_float>(coeffs_rat_imag_(ii,jj)) );
					}
				}
				
				is_rational_coeffs_ = true;

				if(is_hom_vars)
				{
					is_homogenized_ = true;
					hom_variable_ = Integer::Make(0);
				}

			}

			
			
			LinearProduct(VariableGroup const& variables, std::shared_ptr<Node> const& hom_var, Mat<mpfr_complex> const& coeffs, bool is_hom_vars) :
			variables_(variables), num_factors_(coeffs.rows()), hom_variable_(hom_var), is_hom_vars_(is_hom_vars)
			{
				num_variables_ = variables.size();
				
				// Resize coeffs matrix
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				coeffs_dbl_ref.resize(num_factors_, num_variables_+1);
				Mat<mpfr_complex>& coeffs_mpfr_ref = std::get<Mat<mpfr_complex>>(coeffs_);
				coeffs_mpfr_ref.resize(num_factors_, num_variables_+1);
				
				// Resize temporary variable holders
				temp_var_d_.resize(num_variables_ + 1);
				temp_var_d_[num_variables_] = dbl(1);
				temp_var_mp_.resize(num_variables_ + 1);
				temp_var_mp_[num_variables_] = mpfr_complex(1);
				
				coeffs_mpfr_ref = coeffs;
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					for (int jj = 0; jj < num_variables_+1; ++jj)
					{
						coeffs_dbl_ref(ii,jj) = static_cast<dbl>(coeffs_mpfr_ref(ii,jj));
					}
				}
				
				is_rational_coeffs_ = false;
				
				if(is_hom_vars)
				{
					is_homogenized_ = true;
					hom_variable_ = Integer::Make(0);
				}

			}

			
			/**
			 \brief Break off a single linear factor in the product and return as a LinearProduct node.
			 
			 \param index Index of the linear factor, starting at 0
			 \return LinearProduct node contain the single linear.
			 */
			
            std::shared_ptr<LinearProduct> GetLinears(size_t index) const;

			
			
			/**
			 \brief Break off a set of linear factors in the product and return as a LinearProduct node.
			 
			 \param indices std::vector of indices into the factors_ vector
			 \return LinearProduct node contain the linears.
			 */
			
            std::shared_ptr<LinearProduct> GetLinears(std::vector<size_t> indices) const;

			
			
			void SetupVariables(size_t num_factors, VariableGroup const& variables)
			{
				num_variables_ = variables.size();
				
				// Resize coeffs matrix
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				coeffs_dbl_ref.resize(num_factors_, num_variables_+1);
				Mat<mpfr_complex>& coeffs_mpfr_ref = std::get<Mat<mpfr_complex>>(coeffs_);
				coeffs_mpfr_ref.resize(num_factors_, num_variables_+1);
				
				// Resize temporary variable holders
				temp_var_d_.resize(num_variables_ + 1);
				temp_var_d_[num_variables_] = dbl(1);
				temp_var_mp_.resize(num_variables_ + 1);
				temp_var_mp_[num_variables_] = mpfr_complex(1);
				
				hom_variable_ = Integer::Make(1);
			}

			
			
			
			friend class boost::serialization::access;
			
			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<Symbol>(*this);
			}
			
			
			void PrecisionChangeSpecific(unsigned prec) const
			{
				temp_sum_mp_.precision(prec);
				for(auto& v : temp_var_mp_)
				{
					v.precision(prec);
				}
				
				Mat<mpfr_complex> coeffs_ref = std::get<Mat<mpfr_complex>>(coeffs_);
				for(int ii = 0; ii < coeffs_ref.rows(); ++ii)
				{
					for(int jj = 0; jj < coeffs_ref.cols(); ++jj)
					{
						coeffs_ref(ii,jj).precision(prec);
						coeffs_ref(ii,jj).real( static_cast<mpfr_float>(coeffs_rat_real_(ii,jj)) );
						coeffs_ref(ii,jj).imag( static_cast<mpfr_float>(coeffs_rat_imag_(ii,jj)) );
					}
				}
			}
			
		};
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		/**
		 \brief Represents the differential of a single linear.
		 
		 This class represents a linear of differentials.  This is the result of differentiating a single linear node.
		 
		 */
		class DiffLinear : public virtual Symbol, public virtual EnableSharedFromThisVirtual<DiffLinear>
		{
		public:
			BERTINI_DEFAULT_VISITABLE()
			
			virtual ~DiffLinear() = default;
			
			template<typename... Ts> 
			static 
			std::shared_ptr<DiffLinear> Make(Ts&& ...ts){ 
				return std::shared_ptr<DiffLinear>( new DiffLinear(ts...) );
			}

		private:
			
			/**
			 \brief Create a linear differential node. Only a linear, not a product.
			 
			 \param linear The linear that we are differentiating.
			 
			 */
			DiffLinear(std::shared_ptr<LinearProduct> const& linear);
			
		public:
			
			
			
			
			
			
			
			
			
			
			
			
			/**
			 \brief Reset variable values in this node
			 */
			void Reset() const override
			{
				Node::ResetStoredValues();
			};
			
			
			/**
			 Method for printing to output stream
			 */
			void print(std::ostream & target) const override;
			
			
			
			
			/**
			 Return SumOperator whose children are derivatives of children_
			 */
			std::shared_ptr<Node> Differentiate(std::shared_ptr<Variable> const& v = nullptr) const override
			{
				return Integer::Make(0);
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
			
			
			bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
			{
				return true;
			}
			
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
			virtual void precision(unsigned int prec) const override
			{
				auto& val_pair = std::get< std::pair<mpfr_complex,bool> >(current_value_);
				val_pair.first.precision(prec);
				
				this->PrecisionChangeSpecific(prec);
				
			};
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
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
			mpfr_complex FreshEval_mp(std::shared_ptr<Variable> const& diff_variable) const override
			{
				mpfr_complex eval_value;
				
				this->FreshEval_mp(eval_value, diff_variable);
				return eval_value;
			}
			
			
			/**
			 \brief Evaluation of linear product node IN PLACE.  Returns evaluation value.
			 
			 \param evaluation_value The in place variable that stores the evaluation.
			 \param diff_variable Variable that we are differentiating with respect to.  Only for evaluating Jacobians.
			 */
			void FreshEval_mp(mpfr_complex& evaluation_value, std::shared_ptr<Variable> const& diff_variable) const override
			{
				for(int ii = 0; ii < variables_.size(); ++ii)
				{
					if(diff_variable == variables_[ii])
					{
						auto& coeff_ref = std::get<Mat<mpfr_complex>>(coeffs_);
						evaluation_value = coeff_ref(0,ii);
						return;
					}
				}
				
				// If not one of the affine variables
				if(diff_variable == hom_variable_)
				{
					auto& coeff_ref = std::get<Mat<mpfr_complex>>(coeffs_);
					evaluation_value = coeff_ref(0,variables_.size()-1);
					return;
				}
				
				
				// If none of the variables
				evaluation_value = mpfr_complex(0);
				return;
			}
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
		private:
			//			std::vector< std::vector< std::tuple<mpq_rational, dbl, mpfr_complex> > > coeffs_;
			Mat<mpq_rational> coeffs_rat_real_;   ///< Matrix of real rational coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor, with the final column being the constant coefficient.  These rationals can then be downsampled for each data type.
			
			Mat<mpq_rational> coeffs_rat_imag_;   ///< Same as coeffs_rat_real_ but for imaginary portion of the coefficients.
			
			std::tuple< Mat<dbl>, Mat<mpfr_complex> > coeffs_;   ///< Matrix of coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor with the final column being the constant coefficient.  This is a tuple with one matrix for each data type.
			
			VariableGroup variables_; ///< Differentials of variables used in the linear.
			
			std::shared_ptr<Node> hom_variable_; ///< The homogenizing variable for this variable group.  Initially this is set to an Integer = 0.  When we homogenize, this is set to the variable.
			
			
			size_t num_variables_;  ///< The number of variables in each linear.
			
			
			
			
			
			
			
			
			
		private:
			
			DiffLinear() = default;
			
			
			
			friend class boost::serialization::access;
			
			template <typename Archive>
			void serialize(Archive& ar, const unsigned version) {
				ar & boost::serialization::base_object<NaryOperator>(*this);
			}
			
			
			void PrecisionChangeSpecific(unsigned prec) const
			{
				Mat<mpfr_complex> coeffs_ref = std::get<Mat<mpfr_complex>>(coeffs_);
				for(int ii = 0; ii < coeffs_ref.rows(); ++ii)
				{
					for(int jj = 0; jj < coeffs_ref.cols(); ++jj)
					{
						coeffs_ref(ii,jj).precision(prec);
					}
				}
			}
			
		};
		
		
		
		
		
	} // re: namespace node
} // re: namespace bertini




#endif
