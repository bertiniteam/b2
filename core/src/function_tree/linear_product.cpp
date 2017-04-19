//This file is part of Bertini 2.
//
//variable.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//linear_product.cpp is distributed in the hope that it will be useful,
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
// linear_product.cpp:  Declares the class LinearProduct.


/**
 \file linear_product.cpp
 
 \brief Provides the LinearProduct Node class and the DiffLinear node class.
 
 */
#ifndef BERTINI_FUNCTION_TREE_LINPRODUCT_CPP
#define BERTINI_FUNCTION_TREE_LINPRODUCT_CPP

#include "bertini2/function_tree/symbols/linear_product.hpp"


template<typename NumType> using Mat = bertini::Mat<NumType>;

namespace  bertini {
	namespace node{

		
		/////////////////////////////////////////
		//
		//     Linear Product methods
		//
		////////////////////////////////////////
		
		std::shared_ptr<Node> LinearProduct::Differentiate() const
		{
			std::shared_ptr<SumOperator> ret_sum;
			std::vector<size_t> indices;
			// First term of product rule
			std::shared_ptr<MultOperator> temp_mult = std::make_shared<MultOperator>(std::make_shared<DiffLinear>(GetLinears(0)));
			for(int ii = 1; ii < num_factors_; ++ii)
			{
				indices.push_back(ii);
			}
			temp_mult *= GetLinears(indices);
			ret_sum = std::make_shared<SumOperator>(temp_mult, true);
			
			
			// Rest of the factors
			for(int ii = 1; ii < num_factors_; ++ii)
			{
				temp_mult = std::make_shared<MultOperator>(std::make_shared<DiffLinear>(GetLinears(ii)));
				indices.clear();
				for(int jj = 0; jj < num_factors_ && jj != ii; ++jj)
				{
					indices.push_back(jj);
				}
				temp_mult *= GetLinears(indices);
				ret_sum->AddChild(temp_mult,true);
			}
			
			
			
			return ret_sum;
			
			
		}
		
		
		
		
		
		
		int LinearProduct::Degree(std::shared_ptr<Variable> const& v) const
		{
			int deg = 0;
			
			// If v is part of the linear product
			if(std::find(variables_.begin(), variables_.end(), v) != std::end(variables_))
			{
				deg = num_factors_;
			}
			
			return deg;
		}
		
		
		
		
		
		int LinearProduct::Degree(VariableGroup const& vars) const
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
		}

		
		
		
		std::vector<int> LinearProduct::MultiDegree(VariableGroup const& vars) const
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
		
		
		
		
		
		void LinearProduct::Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar)
		{
			// Needed if this is for a homogeneous variable group
			if(is_homogenized_)
			{
				return;
			}
			
			
			// Check if vars is the variable group for this linear product
			bool is_vargroup_same = true; // Check if vars and variables_ are the same
			bool is_v_in_vars = false; // Check if at least one v in variables_ is in vars
			
			
			// Is homvar in vars?
			VariableGroup temp_vars = vars;
			temp_vars.erase(std::remove(temp_vars.begin(), temp_vars.end(), homvar), temp_vars.end());
			
			
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
			
		}
		
		
		
		
		
		bool LinearProduct::IsHomogeneous(std::shared_ptr<Variable> const& v) const
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
		}
		
		
		
		
		
		bool LinearProduct::IsHomogeneous(VariableGroup const& vars) const
		{
			bool is_hom = false;
			
			// Check if vars is the variable group for this linear product
			bool is_vargroup_same = true; // Check if vars and variables_ are the same
			bool is_v_in_vars = false; // Check if at least one v in variables_ is in vars
			for (auto const v : variables_)
			{
				
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
		}
		
		
		
		
		
		
		
		
		
		
		
		std::shared_ptr<LinearProduct> LinearProduct::GetLinears(size_t index) const
		{
			if(is_rational_coeffs_)
			{
				Mat<mpq_rational> temp_real(1,num_variables_+1);
				Mat<mpq_rational> temp_imag(1,num_variables_+1);
				for(int jj = 0; jj < num_variables_+1; ++jj)
				{
					temp_real(0,jj) = coeffs_rat_real_(index,jj);
					temp_imag(0,jj) = coeffs_rat_imag_(index,jj);
				}
				
				LinearProduct temp(variables_, hom_variable_, temp_real, temp_imag, is_hom_vars_);
				return std::make_shared<LinearProduct>(temp);
			}
			else
			{
				Mat<mpfr> temp_mpfr(1,num_variables_+1);
				auto& coeffs_mp_ref = std::get<Mat<mpfr>>(coeffs_);
				for(int jj = 0; jj < num_variables_+1; ++jj)
				{
					temp_mpfr(0,jj) = coeffs_mp_ref(index,jj);
				}
				
				LinearProduct temp(variables_, hom_variable_, temp_mpfr, is_hom_vars_);
				return std::make_shared<LinearProduct>(temp);

			}
			
		}

		
		
		std::shared_ptr<LinearProduct> LinearProduct::GetLinears(std::vector<size_t> indices) const
		{
			if(is_rational_coeffs_)
			{
				Mat<mpq_rational> temp_real(indices.size(),num_variables_+1);
				Mat<mpq_rational> temp_imag(indices.size(),num_variables_+1);
				for(int ii = 0; ii < indices.size(); ++ii)
				{
					for(int jj = 0; jj < num_variables_+1; ++jj)
					{
						temp_real(ii,jj) = coeffs_rat_real_(indices[ii],jj);
						temp_imag(ii,jj) = coeffs_rat_imag_(indices[ii],jj);
					}
				}
				
				LinearProduct temp(variables_, hom_variable_, temp_real, temp_imag, is_hom_vars_);
				return std::make_shared<LinearProduct>(temp);
			}
			else
			{
				Mat<mpfr> temp_mpfr(indices.size(),num_variables_+1);
				auto& coeffs_mp_ref = std::get<Mat<mpfr>>(coeffs_);
				for(int ii = 0; ii < indices.size(); ++ii)
				{
					for(int jj = 0; jj < num_variables_+1; ++jj)
					{
						temp_mpfr(ii,jj) = coeffs_mp_ref(indices[ii],jj);
					}
				}
				
				LinearProduct temp(variables_, hom_variable_, temp_mpfr, is_hom_vars_);
				return std::make_shared<LinearProduct>(temp);
			}
			
		}

		
		
		
		
		
		
		
		
		
		
		
		
		/////////////////////////////////////////
		//
		//     Diff Linear methods
		//
		////////////////////////////////////////
		
		DiffLinear::DiffLinear(std::shared_ptr<LinearProduct> const& linear)
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


	}
}

#endif