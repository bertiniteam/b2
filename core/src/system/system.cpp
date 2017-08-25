//This file is part of Bertini 2.
//
//system.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//system.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with system.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// dani brake, university of wisconsin eau claire


#include "bertini2/system/system.hpp"

template<typename NumType> using Vec = bertini::Vec<NumType>;
template<typename NumType> using Mat = bertini::Mat<NumType>;

BOOST_CLASS_EXPORT(bertini::System)



namespace bertini 
{

	JacobianEvalMethod DefaultJacobianEvalMethod()
	{
		return JacobianEvalMethod::Derivatives;
	}

	bool DefaultAutoSimplify()
	{
		return true;
	}

	void swap(System & a, System & b)
	{
		using std::swap;

		swap(a.ungrouped_variables_,b.ungrouped_variables_);
		swap(a.variable_groups_,b.variable_groups_);
		swap(a.hom_variable_groups_,b.hom_variable_groups_);
		swap(a.homogenizing_variables_,b.homogenizing_variables_);

		swap(a.time_order_of_variable_groups_,b.time_order_of_variable_groups_);

		swap(a.have_path_variable_,b.have_path_variable_);
		swap(a.path_variable_,b.path_variable_);

		swap(a.have_ordering_,b.have_ordering_);
		swap(a.variable_ordering_,b.variable_ordering_);

		swap(a.implicit_parameters_,b.implicit_parameters_);
		swap(a.explicit_parameters_,b.explicit_parameters_);

		swap(a.constant_subfunctions_,b.constant_subfunctions_);
		swap(a.subfunctions_,b.subfunctions_);
		swap(a.functions_,b.functions_);

		swap(a.is_differentiated_,b.is_differentiated_);
		swap(a.jacobian_,b.jacobian_);

		swap(a.space_derivatives_,b.space_derivatives_);
		swap(a.time_derivatives_,b.time_derivatives_);

		swap(a.assume_uniform_precision_,b.assume_uniform_precision_);
		swap(a.jacobian_eval_method_,b.jacobian_eval_method_);

		swap(a.precision_,b.precision_);
		swap(a.is_patched_,b.is_patched_);
		swap(a.patch_,b.patch_);
	}

	// the copy constructor
	System::System(System const& other) : System()
	{
		ungrouped_variables_ = other.ungrouped_variables_;
		variable_groups_  = other.variable_groups_;
		hom_variable_groups_ =  other.hom_variable_groups_;
		homogenizing_variables_ = other.homogenizing_variables_;
		have_path_variable_ = other.have_path_variable_;
		path_variable_ = other.path_variable_;
		implicit_parameters_ = other.implicit_parameters_;
		
		patch_ = other.patch_;
		is_patched_ = other.is_patched_;

		jacobian_ = other.jacobian_;
		space_derivatives_ = other.space_derivatives_;
		time_derivatives_ = other.time_derivatives_;

		is_differentiated_ = other.is_differentiated_;

		assume_uniform_precision_ = other.assume_uniform_precision_;
		jacobian_eval_method_ = other.jacobian_eval_method_;

		time_order_of_variable_groups_ = other.time_order_of_variable_groups_;

		current_variable_values_ = other.current_variable_values_;

		variable_ordering_ = other.variable_ordering_;
		have_ordering_ =  other.have_ordering_;

		precision_ = other.precision_;

		// now to do the members which are not simply copied
		constant_subfunctions_.resize(other.constant_subfunctions_.size());
		for (unsigned ii = 0; ii < constant_subfunctions_.size(); ++ii)
			constant_subfunctions_[ii] = MakeFunction(other.constant_subfunctions_[ii]->entry_node());

		subfunctions_.resize(other.subfunctions_.size());
		for (unsigned ii = 0; ii < subfunctions_.size(); ++ii)
			subfunctions_[ii] = MakeFunction(other.subfunctions_[ii]->entry_node());

		functions_.resize(other.functions_.size());
		for (unsigned ii = 0; ii < functions_.size(); ++ii)
			functions_[ii] = MakeFunction(other.functions_[ii]->entry_node());

		explicit_parameters_.resize(other.explicit_parameters_.size());
		for (unsigned ii = 0; ii < explicit_parameters_.size(); ++ii)
			explicit_parameters_[ii] = MakeFunction(other.explicit_parameters_[ii]->entry_node());
	}

	// the assignment operator
	System& System::operator=(const System & other)
	{
		*this = System(other);
		return *this;
	}


	/////////
	//
	//  getters
	//
	/////////////////


	size_t System::NumFunctions() const
	{
		return functions_.size();
	}


	size_t System::NumVariables() const
	{
		return NumHomVariables() + NumNaturalVariables();
	}

	size_t System::NumNaturalVariables() const
	{
		size_t num_vars = 0;

		for (const auto& iter : variable_groups_)
			num_vars += iter.size();
		for (const auto& iter : hom_variable_groups_)
			num_vars += iter.size();
		num_vars += ungrouped_variables_.size();

		return num_vars;
	}

	size_t System::NumHomVariables() const
	{
		return homogenizing_variables_.size();
	}

	size_t System::NumTotalVariableGroups() const
	{
		return NumVariableGroups() + NumHomVariableGroups();
	}

	size_t System::NumVariableGroups() const
	{
		return variable_groups_.size();
	}

	size_t System::NumUngroupedVariables() const
	{
		return ungrouped_variables_.size();
	}


	size_t System::NumHomVariableGroups() const
	{
		return hom_variable_groups_.size();
	}


	size_t System::NumConstants() const
	{
		return constant_subfunctions_.size();
	}

	size_t System::NumParameters() const
	{
		return explicit_parameters_.size();
	}

	size_t System::NumImplicitParameters() const
	{
		return implicit_parameters_.size();
	}


	size_t System::NumTotalFunctions() const
	{
		return NumFunctions() + NumPatches();
	}


	void System::precision(unsigned new_precision) const
	{
		if (this->assume_uniform_precision_ && new_precision == this->precision_)
			return;

		for (const auto& iter : functions_) {
			iter->precision(new_precision);
		}

		for (const auto& iter : subfunctions_) {
			iter->precision(new_precision);
		}

		for (const auto& iter : explicit_parameters_) {
			iter->precision(new_precision);
		}


		for (const auto& iter :implicit_parameters_) {
			iter->precision(new_precision);
		}

		for (const auto& iter : constant_subfunctions_) {
			iter->precision(new_precision);
		}

		if (is_differentiated_)
		{
			switch (jacobian_eval_method_)
			{
				case JacobianEvalMethod::JacobianNode:
					for (const auto& iter : jacobian_)
						iter->precision(new_precision);
					break;
				case JacobianEvalMethod::Derivatives:
					for (const auto& iter : space_derivatives_)
						iter->precision(new_precision);
					for (const auto& iter : time_derivatives_)
						iter->precision(new_precision);
					break;
				// later, case for straight line program?
			}
			
		}

		if (have_path_variable_)
			path_variable_->precision(new_precision);


		for (const auto& iter : homogenizing_variables_)
			iter->precision(new_precision);

		for (const auto& iter : variable_groups_)
			for (const auto& jter : iter)
				jter->precision(new_precision);

		for (const auto& iter : hom_variable_groups_)
			for (const auto& jter : iter)
				jter->precision(new_precision);

		for (const auto& iter : ungrouped_variables_)
			iter->precision(new_precision);

		using bertini::Precision;
		Precision(std::get<Vec<mpfr> >(current_variable_values_),new_precision);

		if (IsPatched())
			patch_.Precision(new_precision);

		precision_ = new_precision;
	}


	void System::Differentiate() const
	{
		switch (jacobian_eval_method_)
		{
			case JacobianEvalMethod::JacobianNode:
			{
				auto num_functions = NumFunctions();
				jacobian_.resize(num_functions);
				for (int ii = 0; ii < num_functions; ++ii)
					jacobian_[ii] = MakeJacobian(functions_[ii]->Differentiate());
				break;
			}
			case JacobianEvalMethod::Derivatives:
			{
				const auto& vars = this->Variables();
				const auto num_vars = NumVariables();
				const auto num_functions = NumFunctions();

				space_derivatives_.resize(num_functions*num_vars);
				// again, computing these in column major, so staying with one variable at a time.
				for (int jj = 0; jj < num_vars; ++jj)
					for (int ii = 0; ii < num_functions; ++ii)
						space_derivatives_[ii+jj*num_functions] = functions_[ii]->Differentiate(vars[jj]);

				if (HavePathVariable())
				{
					const auto& t = path_variable_;
					time_derivatives_.resize(num_functions);
						for (int ii = 0; ii < num_functions; ++ii)
							time_derivatives_[ii] = functions_[ii]->Differentiate(t);
				}
				break;
			}
		}
		is_differentiated_ = true;

		if (auto_simplify_)
		{
			this->SimplifyDerivatives();
		}
	}





	void System::Homogenize()
	{

		// first some checks to make sure the system is compatible with the act of homogenization
		//
		//  a system must be:
		//    * homogeneous already with respect to the homogeneous variable groups
		//    * polynomial
		//    * not partially homogenized, in the sense that some groups have been homogenized, and others haven't
		//    
		//
		for (const auto& curr_function : functions_)
		{	
			for (const auto& curr_var_gp : hom_variable_groups_)
			{
				if (!curr_function->IsHomogeneous(curr_var_gp))
					throw std::runtime_error("inhomogeneous function, with homogeneous variable group");
			}
		}

		if (!IsPolynomial())
			throw std::runtime_error("trying to homogenize a non-polynomial system.");

		bool already_had_homvars = NumHomVariables()!=0;
		
		if (already_had_homvars && NumHomVariables()!=NumVariableGroups())
			throw std::runtime_error("size mismatch on number of homogenizing variables and number of variable groups");

		if (!already_had_homvars)
		{
			homogenizing_variables_.resize(NumVariableGroups());
		}


		auto group_counter = 0;
		for (auto curr_var_gp = variable_groups_.begin(); curr_var_gp!=variable_groups_.end(); curr_var_gp++)
		{
			std::stringstream converter;
			converter << "HOM_VAR_" << group_counter;

			if (already_had_homvars){
				Var hom_var = homogenizing_variables_[group_counter];
				VariableGroup temp_group = *curr_var_gp;
				temp_group.push_front(hom_var);
				for (const auto& curr_function : functions_)
					curr_function->Homogenize(temp_group, hom_var);
			}
			else
			{
				Var hom_var = MakeVariable(converter.str());
				homogenizing_variables_[group_counter] = hom_var;
				for (const auto& curr_function : functions_)
					curr_function->Homogenize(*curr_var_gp, hom_var);
			}

			group_counter++;
		}

		is_differentiated_ = false;
		have_ordering_ = false;

		#ifndef BERTINI_DISABLE_ASSERTS
		assert(homogenizing_variables_.size() == variable_groups_.size());
		#endif
	}




	bool System::IsHomogeneous() const
	{
		bool have_homvars = NumHomVariables()!=0;

		if (NumHomVariables()!=NumVariableGroups())
			return false;

		for (const auto& iter : functions_)
		{
			auto counter = 0;
			for (const auto& vars : variable_groups_)
			{
				auto tempvars = vars;
				if (have_homvars)
					tempvars.push_front(homogenizing_variables_[counter]);
				counter++;

				if (!iter->IsHomogeneous(tempvars))
					return false;
			}

			for (const auto& vars : hom_variable_groups_)
				if (!iter->IsHomogeneous(vars))
					return false;

			if (NumUngroupedVariables()>0)
				if (!iter->IsHomogeneous(ungrouped_variables_))
					return false;
		}
		return true;
	}


	bool System::IsPolynomial() const
	{
		bool have_homvars = NumHomVariables()!=0;
		if (have_homvars && NumHomVariables()!=NumVariableGroups())
			throw std::runtime_error("trying to check polynomiality on a partially-formed system.  mismatch between number of homogenizing variables, and number of variable groups");


		for (const auto& iter : functions_)
		{
			auto counter = 0;
			for (const auto& vars : variable_groups_)
			{
				auto tempvars = vars;
				if (have_homvars)
					tempvars.push_front(homogenizing_variables_[counter]);

				counter++;

				if (!iter->IsPolynomial(tempvars))
					return false;

			}
			for (const auto& vars : hom_variable_groups_)
				if (!iter->IsPolynomial(vars))
					return false;

		}
		return true;
	}








	////////////////////
	//
	//  Adders
	//
	//////////////////////





	void System::AddVariableGroup(VariableGroup const& v)
	{
		variable_groups_.push_back(v);
		is_differentiated_ = false;
		have_ordering_ = false;
		is_patched_ = false;
		time_order_of_variable_groups_.push_back( VariableGroupType::Affine);
	}




	void System::AddHomVariableGroup(VariableGroup const& v)
	{
		hom_variable_groups_.push_back(v);
		is_differentiated_ = false;
		have_ordering_ = false;
		is_patched_ = false;
		time_order_of_variable_groups_.push_back( VariableGroupType::Homogeneous);
	}





	void System::AddUngroupedVariable(Var const& v)
	{
		ungrouped_variables_.push_back(v);
		is_differentiated_ = false;
		have_ordering_ = false;
		is_patched_ = false;
		time_order_of_variable_groups_.push_back( VariableGroupType::Ungrouped);
	}




	void System::AddUngroupedVariables(VariableGroup const& v)
	{
		ungrouped_variables_.insert( ungrouped_variables_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
		have_ordering_ = false;
		is_patched_ = false;
		for (const auto& iter : v)
			time_order_of_variable_groups_.push_back( VariableGroupType::Ungrouped);
	}



 
	void System::AddImplicitParameter(Var const& v)
	{
		implicit_parameters_.push_back(v);
		is_differentiated_ = false;
	}




	void System::AddImplicitParameters(VariableGroup const& v)
	{
		implicit_parameters_.insert( implicit_parameters_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
	}









	void System::AddParameter(Fn const& F)
	{
		explicit_parameters_.push_back(F);
		is_differentiated_ = false;
	}



	void System::AddParameters(std::vector<Fn> const& v)
	{
		explicit_parameters_.insert( explicit_parameters_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
	}





	void System::AddSubfunction(Fn const& F)
	{
		subfunctions_.push_back(F);
		is_differentiated_ = false;
	}



	void System::AddSubfunctions(std::vector<Fn> const& v)
	{
		subfunctions_.insert( subfunctions_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
	}





	void System::AddFunction(Fn const& F)
	{
		functions_.push_back(F);
		is_differentiated_ = false;
	}



	void System::AddFunction(Nd const& N)
	{
		Fn F = MakeFunction(N);
		functions_.push_back(F);
		is_differentiated_ = false;
	}



	void System::AddFunctions(std::vector<Fn> const& v)
	{
		functions_.insert( functions_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
	}






	void System::AddConstant(Fn const& F)
	{
		constant_subfunctions_.push_back(F);
		is_differentiated_ = false;
	}


	void System::AddConstants(std::vector<Fn> const& v)
	{
		constant_subfunctions_.insert( constant_subfunctions_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
	}






	void System::AddPathVariable(Var const& v)
	{
		path_variable_ = v;
		is_differentiated_ = false;
		have_path_variable_ = true;
	}



	bool System::HavePathVariable() const
	{
		return have_path_variable_;
	}

	










	/////////////////
	//
	// Variable ordering functions
	//
	///////////////////

	VariableGroup System::VariableOrdering() const
	{
		bool have_homvars = NumHomVariables()>0;
		if (have_homvars && NumHomVariables() != NumVariableGroups())
			throw std::runtime_error("mismatch between number of homogenizing variables, and number of affine variables groups.  unable to form variable vector in FIFO ordering.  you probably need to homogenize the system.");

		

		VariableGroup constructed_ordering;

		unsigned affine_group_counter = 0, ungrouped_variable_counter = 0, hom_group_counter = 0;
		for (auto group_type : time_order_of_variable_groups_)
		{
			switch (group_type){
				case VariableGroupType::Affine:
				{
					if (have_homvars)
						constructed_ordering.push_back(homogenizing_variables_[affine_group_counter]);

					constructed_ordering.insert(constructed_ordering.end(), 
					                variable_groups_[affine_group_counter].begin(), 
					                variable_groups_[affine_group_counter].end());
					affine_group_counter++;
					break;
				}
				case VariableGroupType::Homogeneous:
				{
					constructed_ordering.insert(constructed_ordering.end(), 
					                hom_variable_groups_[hom_group_counter].begin(), 
					                hom_variable_groups_[hom_group_counter].end());
					hom_group_counter++;
					break;
				}
				case VariableGroupType::Ungrouped:
				{
					constructed_ordering.push_back(ungrouped_variables_[ungrouped_variable_counter]);
					ungrouped_variable_counter++;
					break;
				}
				default:
				{	
					throw std::runtime_error("unacceptable VariableGroupType in FIFOVariableOrdering");
				}
			}
		}
		
		#ifndef BERTINI_DISABLE_ASSERTS
		assert(constructed_ordering.size()==NumVariables() && "resulting constructed ordering has differing size from the number of variables in the problem.");
		#endif

		return constructed_ordering;	
	}


	

	//private
	void System::ConstructOrdering() const
	{
		variable_ordering_ = VariableOrdering();
		have_ordering_ = true;
	}



	const VariableGroup& System::Variables() const
	{
		if (!have_ordering_){
			ConstructOrdering();
		}

		return variable_ordering_;
	}
		

	void System::CopyVariableStructure(System const& other)
	{
		this->ClearVariables();

		time_order_of_variable_groups_ = other.time_order_of_variable_groups_;

		this->ungrouped_variables_ = other.ungrouped_variables_;
		this->variable_groups_ = other.variable_groups_;
		hom_variable_groups_ = other.hom_variable_groups_;
		homogenizing_variables_ = other.homogenizing_variables_;

		path_variable_ = other.path_variable_;
		have_path_variable_ = other.have_path_variable_;

		variable_ordering_ = other.variable_ordering_; 
		have_ordering_ = other.have_ordering_;
	}





	std::vector<unsigned> System::VariableGroupSizesFIFO() const
	{
		bool have_homvars = NumHomVariables()>0;
		if (have_homvars && NumHomVariables() != NumVariableGroups())
			throw std::runtime_error("mismatch between number of homogenizing variables, and number of affine variables groups.  unable to form variable vector in FIFO ordering.  you probably need to homogenize the system.");

		std::vector<unsigned> s;

		unsigned hom_group_counter(0), affine_group_counter(0), patch_counter(0);
		for (auto curr_grouptype : time_order_of_variable_groups_)
		{
			if (curr_grouptype==VariableGroupType::Homogeneous)
				s.push_back(hom_variable_groups_[hom_group_counter++].size());
			else if (curr_grouptype==VariableGroupType::Affine)
				s.push_back(variable_groups_[affine_group_counter++].size() + static_cast<int>(have_homvars));
		}
		return s;
	}


	/////////////////
	//
	// Patching functions
	//
	///////////////////
	


	void System::AutoPatchFIFO()
	{
		if (!IsHomogeneous())
			throw std::runtime_error("requesting to AutoPatch a system which is not homogenized.  Homogenize it first.");
		
		patch_ = Patch(VariableGroupSizesFIFO());

		is_patched_ = true;
	}



	void System::CopyPatches(System const& other)
	{
		if (!other.IsPatched())
			throw std::runtime_error("trying to copy patch from unpatched other system.  may only copy patch from a system which is already patched.");

		this->patch_ = other.patch_;
		is_patched_ = true;
	}


			



    /////////////////
	//
	// Functions involving the coefficients and degrees of functions in the systems.
	//
	///////////////////

	template <typename NumT>
	typename Eigen::NumTraits<NumT>::Real System::CoefficientBound(unsigned num_evaluations) const
	{
		static_assert(Eigen::NumTraits<NumT>::IsComplex,"NumT must be a complex type");
		
		using RT = typename Eigen::NumTraits<NumT>::Real;
		using CT = NumT;

		RT bound(0);

		for (unsigned ii=0; ii < num_evaluations; ii++)
		{	
			Vec<CT> randy = RandomOfUnits<CT>(NumVariables());
			Vec<CT> f_vals;
			if (HavePathVariable())
				f_vals = Eval(randy, RandomUnit<CT>());
			else
				f_vals = Eval(randy);
			
			Mat<CT> dh_dx = Jacobian<CT>();
			
			bound = max(f_vals.array().abs().maxCoeff(),
						 dh_dx.array().abs().maxCoeff(), bound);
		}
		return bound;
	}

	template double System::CoefficientBound<dbl>(unsigned) const;
	template mpfr_float System::CoefficientBound<mpfr>(unsigned) const;

    int System::DegreeBound() const
    {
    	auto degs = Degrees(Variables());
    	return *std::max_element(degs.begin(), degs.end());
    }


	std::vector<int> System::Degrees() const
	{
		std::vector<int> degs;
		for (const auto& iter : functions_)
			degs.push_back(iter->Degree());
		return degs;
	}


	std::vector<int> System::Degrees(VariableGroup const& vars) const
	{
		std::vector<int> degs;
		for (const auto& iter : functions_)
			degs.push_back(iter->Degree(vars));
		return degs;
		}


	void System::ReorderFunctionsByDegreeDecreasing()
	{
		auto degs = Degrees(Variables());

		// now we sort a vector of the indexing numbers by the degrees contained in degs.
		std::vector<size_t> indices(degs.size());
		//http://en.cppreference.com/w/cpp/algorithm/iota
		//http://www.cplusplus.com/doc/tutorial/typecasting/
		std::iota(begin(indices), end(indices), static_cast<size_t>(0));
		std::sort( begin(indices), end(indices), [&](size_t a, size_t b) { return degs[a] > degs[b]; } );
		


		// finally, we re-order the functions based on the indices we just computed
		std::vector<std::shared_ptr<node::Function> > re_ordered_functions(degs.size());
		size_t ind = 0;
		for (auto iter : indices)
		{
			re_ordered_functions[ind] = functions_[iter];
			ind++;
		}

		swap(functions_, re_ordered_functions);
	}



	void System::ReorderFunctionsByDegreeIncreasing()
	{
		auto degs = Degrees(Variables());

		// now we sort a vector of the indexing numbers by the degrees contained in degs.
		std::vector<size_t> indices(degs.size());
		//http://en.cppreference.com/w/cpp/algorithm/iota
		//http://www.cplusplus.com/doc/tutorial/typecasting/
		std::iota(begin(indices), end(indices), static_cast<size_t>(0));
		std::sort( begin(indices), end(indices), [&](size_t a, size_t b) { return degs[a] < degs[b]; } );
		


		// finally, we re-order the functions based on the indices we just computed
		std::vector<std::shared_ptr<node::Function> > re_ordered_functions(degs.size());
		size_t ind = 0;
		for (auto iter : indices)
		{
			re_ordered_functions[ind] = functions_[iter];
			ind++;
		}

		swap(functions_, re_ordered_functions);
	}











	/////////////////
	//
	// Clearing functions
	//
	///////////////////

	void System::ClearVariables()
	{
		ungrouped_variables_.clear();
		variable_groups_.clear();
		hom_variable_groups_.clear();
		homogenizing_variables_.clear();

		path_variable_.reset();
		have_path_variable_ = false;
	}



	void System::SimplifyFunctions()
	{
		using bertini::Simplify;
		for (auto& iter : this->functions_)
			Simplify(iter);
	}



	void System::SimplifyDerivatives() const
	{
		using bertini::Simplify;

		auto num_vars = this->NumVariables();
		std::vector<dbl> old_vals(num_vars);  dbl old_path_var_val;

		auto vars = this->Variables();
		for (unsigned ii=0; ii<num_vars; ++ii)
		{
			old_vals[ii] = vars[ii]->Eval<dbl>();
			vars[ii]->SetToRandUnit<dbl>();
		}

		if (HavePathVariable())
		{
			old_path_var_val = path_variable_->Eval<dbl>();
			path_variable_->SetToRandUnit<dbl>();
		}


		for (const auto& n : jacobian_)
			n->Reset();
		for (const auto& n : space_derivatives_)
			n->Reset();
		for (const auto& n : time_derivatives_)
			n->Reset();


		switch (jacobian_eval_method_)
		{
			case JacobianEvalMethod::JacobianNode:
				for (auto& iter : this->jacobian_)
					Simplify(iter);
				break;
			case JacobianEvalMethod::Derivatives:
				for (auto& iter : this->space_derivatives_)
					Simplify(iter);
				for (auto& iter : this->time_derivatives_)
					Simplify(iter);
				break;

		}
		
		for (unsigned ii=0; ii<num_vars; ++ii)
			vars[ii]->set_current_value<dbl>(old_vals[ii]);
		if (HavePathVariable())
			path_variable_->set_current_value(old_path_var_val);


		for (const auto& n : jacobian_)
			n->Reset();
		for (const auto& n : space_derivatives_)
			n->Reset();
		for (const auto& n : time_derivatives_)
			n->Reset();

	}



	void System::Simplify()
	{
		SimplifyFunctions();
		SimplifyDerivatives();
	}










	
	//////////////////
	//
	//  output operators
	//
	////////////////////

	std::ostream& operator<<(std::ostream& out, const bertini::System & s)
	{


		out << s.NumVariableGroups() << " variable groups, containing these variables:\n";
		auto counter = 0;
		for (const auto& iter : s.variable_groups_)
		{
			out << "group " << counter << ": "<< "\n";
			for (auto jter : iter)
				out << *jter << " ";


			out << "\n";
			counter++;
		}

		out << s.NumHomVariables() << " homogenizing variables:\n";
		for (const auto& iter : s.homogenizing_variables_)
			out << (*iter) << " ";
		out << "\n";


		out << s.NumFunctions() << " functions:\n";
		for (const auto& iter : s.functions_) 
			out << (iter)->name() << " = " << *iter << "\n";
		out << "\n";


		if (s.NumParameters()) {
			out << s.NumParameters() << " explicit parameters:\n";
			for (const auto& iter : s.explicit_parameters_)
				out << (iter)->name() << " = " << *iter << "\n";
			out << "\n";
		}


		if (s.NumConstants()) {
			out << s.NumConstants() << " constants:\n";
			for (const auto& iter : s.constant_subfunctions_)
				out << (iter)->name() << " = " << *iter << "\n";
			out << "\n";
		}

		if (s.path_variable_)
			out << "path variable defined.  named " << s.path_variable_->name() << "\n";
		else 
			out << "no path variable defined\n";

		if (s.is_differentiated_)
		{
			out << "system is differentiated; jacobian:\n";
			switch (s.jacobian_eval_method_)
			{
			case JacobianEvalMethod::JacobianNode:
				for (const auto& iter : s.jacobian_)
					out << (iter)->name() << " = " << *iter << "\n";
				break;
			case JacobianEvalMethod::Derivatives:
				for (int jj = 0; jj < s.NumVariables(); ++jj)
					for (int ii = 0; ii < s.NumFunctions(); ++ii)
					{
						const auto& d = s.space_derivatives_[ii+jj*s.NumFunctions()];
						out << "jac_space_der(" << ii << "," << jj << ") = " << d << "\n";
					}

				if (s.HavePathVariable())
					for (int ii = 0; ii < s.NumFunctions(); ++ii)
					{
						const auto& d = s.time_derivatives_[ii];
						out << "jac_time_der(" << ii << ") = " << d << "\n";
					}
				break;
			}
			out << "\n";
		}
		else{
			out << "system not differentiated\n";
		}

		if (s.IsPatched())
		{
			out << s.patch_;
		}
		else{
			out << "system not patched\n";
		}

		return out;
	}









	/////////////////
	//
	// Arithemetic operators
	//
	///////////////////


	System& System::operator+=(System const& rhs)
	{
		if (this->NumFunctions()!=rhs.NumFunctions())
			throw std::runtime_error("cannot add two Systems with differing numbers of functions");

		if (this->NumVariables()!=rhs.NumVariables())
			throw std::runtime_error("cannot add two Systems with differing numbers of variables");

		if (this->NumHomVariables()!=rhs.NumHomVariables())
			throw std::runtime_error("cannot add two Systems with differing numbers of homogenizing variables");

		if (this->NumTotalVariableGroups()!=rhs.NumTotalVariableGroups())
			throw std::runtime_error("cannot add two Systems with differing total numbers of variable groups");


		//
		//  deal with the patches
		//
		if (!this->IsPatched() && rhs.IsPatched())
			CopyPatches(rhs);

		// the condition (this->IsPatched() && !rhs.IsPatched()) is ok.  nothing to do.
		// the condition (!this->IsPatched() && !rhs.IsPatched()) is ok.  nothing to do.
		else if (this->IsPatched() && rhs.IsPatched())
			if (this->patch_ != rhs.patch_)
				throw std::runtime_error("System+=System cannot combine two patched systems whose patches differ.");

		for (auto iter=functions_.begin(); iter!=functions_.end(); iter++)
			(*iter)->SetRoot( (*(rhs.functions_.begin()+(iter-functions_.begin())))->entry_node() + (*iter)->entry_node());

		return *this;
	}

	const System operator+(System lhs, System const& rhs)
	{
		return lhs+=rhs;
	}


	System& System::operator*=(std::shared_ptr<node::Node> const& N)
	{
		for (auto iter=functions_.begin(); iter!=functions_.end(); iter++)
		{
			(*iter)->SetRoot( N * (*iter)->entry_node());
		}
		return *this;
	}


	const System operator*(System s, std::shared_ptr<node::Node> const&  N)
	{
		return s*=N;
	}


	const System operator*(std::shared_ptr<node::Node> const&  N, System const& s)
	{
		return s*N;
	}


	
	System Concatenate(System sys1, System const& sys2)
	{
		// first we will deal with the variable structure
		if (sys1.NumVariables()!=sys2.NumVariables())
			throw std::runtime_error("concatenating systems with differing numbers of variables");

		if (sys1.VariableOrdering() != sys2.VariableOrdering())
			throw std::runtime_error("concatenating systems with differing variable orderings");

		if (sys1.IsPatched() && sys2.IsPatched())
			if (sys1.GetPatch()!=sys2.GetPatch())
				throw std::runtime_error("concatenating systems with incompatible patches");

		if (sys2.IsPatched() && !sys1.IsPatched())
			sys1.CopyPatches(sys1);
		// the other cases are automatically covered.  sys1 already patched, or neither patched.

		for (unsigned ii(0); ii<sys2.NumFunctions(); ++ii)
			sys1.AddFunction(sys2.Function(ii));

		return sys1;
	}


	System Clone(System const& sys)
	{

//////////////////  attempt 1.  generates a npos == null problem of some sort.  i couldn't figure it out.


		// namespace io = boost::iostreams;
		// using buffer_type = std::vector<char>;
		// buffer_type buffer;

		// io::stream<io::back_insert_device<buffer_type> > output_stream(buffer);
		// boost::archive::binary_oarchive oa(output_stream);

		// oa << sys;
		// output_stream.flush();

		

		// io::basic_array_source<char> source(&buffer[0],buffer.size());
		// io::stream<io::basic_array_source <char> > input_stream(source);
		// boost::archive::binary_iarchive ia(input_stream);

		// System sys_clone;
		// ia >> sys_clone;

		// return sys_clone;



///////////////////////  attempt2  generates crashes.  :(
		// std::string serial_str;
		// {
		// 	boost::iostreams::back_insert_device<std::string> inserter(serial_str);
		// 	boost::iostreams::stream<boost::iostreams::back_insert_device<std::string> > s(inserter);
		// 	boost::archive::binary_oarchive oa(s);

		// 	oa << sys;

		// 	// don't forget to flush the stream to finish writing into the buffer
		// 	s.flush();
		// }
		
		// boost::iostreams::basic_array_source<char> device(serial_str.data(), serial_str.size());
		// boost::iostreams::stream<boost::iostreams::basic_array_source<char> > t(device);
		// boost::archive::binary_iarchive ia(t);
		// System sys_clone;
		// ia >> sys_clone;




///////////////////// attempt3.  works.  why the others generate problems with the binary archive baffles me.

		std::stringstream ss;
		{
			boost::archive::text_oarchive oa(ss);
			oa << sys;
		}

		System sys_clone;
		{
			boost::archive::text_iarchive ia(ss);
			ia >> sys_clone;
		}

		return sys_clone;
	}

	void Simplify(System & sys)
	{
		sys.Simplify();
	}

}
