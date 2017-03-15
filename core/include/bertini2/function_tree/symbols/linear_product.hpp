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
		class LinearProduct : public virtual Symbol
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
				factors_.resize(num_factors_);
				
				// Resize coeffs matrix
				coeffs_rat_real_.resize(num_factors_, num_variables_ + 1);
				coeffs_rat_imag_.resize(num_factors_, num_variables_ + 1);
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				coeffs_dbl_ref.resize(num_factors_, num_variables_ + 1);
				Mat<mpfr>& coeffs_mpfr_ref = std::get<Mat<mpfr>>(coeffs_);
				coeffs_mpfr_ref.resize(num_factors_, num_variables_ + 1);
				
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					for (int jj = 0; jj < num_variables_ + 1; ++jj)
					{
						// Generate random constants as mpq_rationals.  Then downsample to mpfr and dbl.
						// TODO: RandomRat() does not generate random numbers.  Same each run.
						coeffs_rat_real_(ii,jj) = RandomRat();
						coeffs_rat_imag_(ii,jj) = RandomRat();
						coeffs_dbl_ref(ii,jj).real( static_cast<double>(coeffs_rat_real_(ii,jj)) );
						coeffs_dbl_ref(ii,jj).imag( static_cast<double>(coeffs_rat_imag_(ii,jj)) );
						coeffs_mpfr_ref(ii,jj).real( static_cast<mpfr_float>(coeffs_rat_real_(ii,jj)) );
						coeffs_mpfr_ref(ii,jj).imag( static_cast<mpfr_float>(coeffs_rat_imag_(ii,jj)) );
						
						// Create i^th factor as sumoperator and add term
						if(jj == 0)
						{
							factors_[ii] =
								std::make_shared<SumOperator>(std::make_shared<MultOperator>(bertini::MakeRational(coeffs_rat_real_(ii,jj), coeffs_rat_imag_(ii,jj)), variables_[jj]), true);
						}
						else if(jj < num_variables_)
						{
							factors_[ii]->AddChild(std::make_shared<MultOperator>(bertini::MakeRational(coeffs_rat_real_(ii,jj), coeffs_rat_imag_(ii,jj)), variables_[jj]), true);
						}
						else
						{
							factors_[ii]->AddChild(bertini::MakeRational( coeffs_rat_real_(ii,jj), coeffs_rat_imag_(ii,jj) ), true);
						}
					} //re: variable loop
				}//re: factor loop
			}
			

			
			/**
			 /brief Constructor for a linear product node that passes in random coefficients.
			 
			 \param variables A deque of variable nodes that are used in each linear factor.  This does not have to be an actual system variable group.
			 \param num_factors The number of linear factors in the product.
			 \param coeffs_real A matrix of dbl data types for all the coefficients in the linear factors.  Matrix must be of size num_factors x (num_variables + 1) with constant coefficient in last column.
			 \param coeffs_mpfr A matrix of mpfr data types for all the coefficients in the linear factors.  Matrix must be of size num_factors x (num_variables + 1) with constant coefficient in last column.
			 
			 */
			LinearProduct(VariableGroup variables, int num_factors, Mat<dbl>& coeffs_dbl, Mat<mpfr>& coeffs_mpfr)
					: variables_(variables), num_factors_(num_factors)
			{
				num_variables_ = variables.size();
				factors_.resize(num_factors_);
				
				// Resize coeffs matrix
				Mat<dbl>& coeffs_dbl_ref = std::get<Mat<dbl>>(coeffs_);
				coeffs_dbl_ref.resize(num_factors_, num_variables_);
				Mat<mpfr>& coeffs_mpfr_ref = std::get<Mat<mpfr>>(coeffs_);
				coeffs_mpfr_ref.resize(num_factors_, num_variables_);
				
				
				// Set the coefficient matrices with input matrices.
				coeffs_dbl_ref = coeffs_dbl;
				coeffs_mpfr_ref = coeffs_mpfr;
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					for (int jj = 0; jj < num_variables_ + 1; ++jj)
					{
						// Create i^th factor as sumoperator and add term
						if(jj == 0)
						{
							factors_[ii] =
							std::make_shared<SumOperator>(std::make_shared<MultOperator>(bertini::MakeFloat(coeffs_mpfr_ref(ii,jj)), variables_[jj]), true);
						}
						else if(jj < num_variables_)
						{
							factors_[ii]->AddChild(std::make_shared<MultOperator>(bertini::MakeFloat(coeffs_mpfr_ref(ii,jj)), variables_[jj]), true);
						}
						else
						{
							factors_[ii]->AddChild(bertini::MakeFloat(coeffs_mpfr_ref(ii,jj)) , true);
						}
					} //re: variable loop
				}//re: factor loop
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
				
				for (auto f = factors_.begin(); f != factors_.end(); ++f)
				{
					auto factor_deg = (*f)->Degree(v);
					
					deg += factor_deg;
				}
				
//				// If v is part of the linear product
//				if(std::find(variables_.begin(), variables_.end(), v) != std::end(variables_))
//				{
//					deg = num_factors_;
//				}
//				
				
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
				
				for (auto f : factors_)
				{
					auto factor_deg = f->Degree(vars);
					
					deg += factor_deg;
				}

//				
//				for (auto v = vars.begin(); v != vars.end(); ++v)
//				{
//					// if v is a part of the linear product
//					if(std::find(variables_.begin(), variables_.end(), *v) != std::end(variables_))
//					{
//						deg = num_factors_;
//						break;
//					}
//				}
				
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
				
				for (f : factors_)
				{
					auto factor_deg = f->MultiDegree(vars);
					
					for (auto iter = factor_deg.begin(); iter != factor_deg.end(); ++iter)
					{
						*(degs.begin()+(iter-factor_deg.begin())) += *iter;
					}
				}
				
				
//				for (auto v = vars.begin(); v != vars.end(); ++v)
//				{
//					// If v is part of the linear product
//					if(std::find(variables_.begin(), variables_.end(), *v) != std::end(variables_))
//					{
//						*(degs.begin()+(v-vars.begin())) = num_factors_;
//					}
//				}
				
				
				return degs;
			}
			
			
			
			
			
			/**
			 \brief Homogenize a sum, with respect to a variable group, and using a homogenizing variable.
			 
			 \param vars Variable group to homogenize with respect to.
			 \param homvar Homogenization variable.
			 */
			void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
			{
				
			};
			
			bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override{return false;};
			
			/**
			 Check for homogeneity, with respect to a variable group.
			 */
			bool IsHomogeneous(VariableGroup const& vars) const override {return false;};
			
			
			
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
				
				
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					factors_[ii]->EvalInPlace(temp_d_, diff_variable);
					
					// Multiply factors together
					evaluation_value *= temp_d_;
				}
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
				
				for (int ii = 0; ii < num_factors_; ++ii)
				{
					factors_[ii]->EvalInPlace(temp_mp_, diff_variable);
					
					// Multiply factors together
					evaluation_value *= temp_mp_;
				}
			}

			
			
			
			
			
		private:
//			std::vector< std::vector< std::tuple<mpq_rational, dbl, mpfr> > > coeffs_;
			Mat<mpq_rational> coeffs_rat_real_;   ///< Matrix of real rational coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor, with the final column being the constant coefficient.  These rationals can then be downsampled for each data type.
			Mat<mpq_rational> coeffs_rat_imag_;   ///< Same as coeffs_rat_real_ but for imaginary portion of the coefficients.
			std::tuple< Mat<dbl>, Mat<mpfr> > coeffs_;   ///< Matrix of coefficients that define the linear product.  Each row corresponds to a factor in the product, columns correspond to the terms in each factor with the final column being the constant coefficient.  This is a tuple with one matrix for each data type.
			VariableGroup variables_; ///< Variables to be used in each linear factor.  Does not have to correspond directly to a variable group from the system.
			
			std::vector< std::shared_ptr<SumOperator> > factors_;  ///< All the various linear factors in the product.
			
			size_t num_factors_;  ///< The number of factors in the linear product.
			size_t num_variables_;  ///< The number of variables in each linear.
			
			
			mutable mpfr temp_mp_;
			mutable dbl temp_d_;
			mutable mpfr temp_sum_mp_;
			mutable dbl temp_sum_d_;
			
			
			
			
			
			
			
		private:
			
			LinearProduct() = default;
			
			
			
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
