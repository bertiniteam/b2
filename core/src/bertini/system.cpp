//This file is part of Bertini 2.0.
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

//  system.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015


#include "system.hpp"

BOOST_CLASS_EXPORT(bertini::System)

template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

namespace bertini 
{



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
		size_t num_vars = 0;

		num_vars+=homogenizing_variables_.size();
		for (auto iter : variable_groups_)
			num_vars+= iter.size();
		for (auto iter : hom_variable_groups_)
			num_vars+= iter.size();
		num_vars+=ungrouped_variables_.size();

		return num_vars;
	}

	size_t System::NumNaturalVariables() const
	{
		size_t num_vars = 0;

		for (auto iter : variable_groups_)
			num_vars+= iter.size();
		for (auto iter : hom_variable_groups_)
			num_vars+= iter.size();
		num_vars+=ungrouped_variables_.size();

		return num_vars;
	}

	size_t System::NumHomVariables() const
	{
		return homogenizing_variables_.size();
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


	/**
	 Get the number of implicit parameters in this system
	 */
	size_t System::NumImplicitParameters() const
	{
		return implicit_parameters_.size();
	}




	void System::precision(unsigned new_precision)
	{
		for (auto iter : functions_) {
			iter->precision(new_precision);
		}

		for (auto iter : subfunctions_) {
			iter->precision(new_precision);
		}

		for (auto iter : explicit_parameters_) {
			iter->precision(new_precision);
		}


		for (auto iter :implicit_parameters_) {
			iter->precision(new_precision);
		}

		for (auto iter : constant_subfunctions_) {
			iter->precision(new_precision);
		}

		if (is_differentiated_)
			for (auto iter : jacobian_)
				iter->precision(new_precision);

		if (have_path_variable_)
			path_variable_->precision(new_precision);


		for (auto iter : homogenizing_variables_)
			iter->precision(new_precision);

		for (auto iter : variable_groups_)
			for (auto jter : iter)
				jter->precision(new_precision);

		for (auto iter : hom_variable_groups_)
			for (auto jter : iter)
				jter->precision(new_precision);

		for (auto iter : ungrouped_variables_)
			iter->precision(new_precision);


		precision_ = new_precision;
	}




	template<typename T>
	Vec<T> System::Eval(const Vec<T> & variable_values)
	{

		if (variable_values.size()!=NumVariables())
			throw std::runtime_error("trying to evaluate system, but number of variables doesn't match.");
		if (have_path_variable_)
			throw std::runtime_error("not using a time value for evaluation of system, but path variable IS defined.");


		// this function call traverses the entire tree, resetting everything.
		//
		// TODO: it has the unfortunate side effect of resetting constant functions, too.
		//
		// we need to work to correct this.
		for (auto iter : functions_) {
			iter->Reset();
		}

		SetVariables(variable_values);


		Vec<T> value(NumFunctions()); // create vector with correct number of entries.

		{ // for scoping of the counter.
			auto counter = 0;
			for (auto iter=functions_.begin(); iter!=functions_.end(); iter++, counter++) {
				value(counter) = (*iter)->Eval<T>();
			}
		}

		return value;
	}

	// these two lines are explicit instantiations of the template above.  template definitions separate from declarations cause linking problems.  
	// see
	// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
	template Vec<dbl> System::Eval(const Vec<dbl> & variable_values);
	template Vec<mpfr> System::Eval(const Vec<mpfr> & variable_values);


	template<typename T>
	Vec<T> System::Eval(const Vec<T> & variable_values, const T & path_variable_value)
	{

		if (variable_values.size()!=NumVariables())
			throw std::runtime_error("trying to evaluate system, but number of variables doesn't match.");
		if (!have_path_variable_)
			throw std::runtime_error("trying to use a time value for evaluation of system, but no path variable defined.");


		// this function call traverses the entire tree, resetting everything.
		//
		// TODO: it has the unfortunate side effect of resetting constant functions, too.
		//
		// we need to work to correct this.
		for (auto iter : functions_) {
			iter->Reset();
		}

		SetVariables(variable_values);
		SetPathVariable(path_variable_value);


		Vec<T> value(NumFunctions()); // create vector with correct number of entries.

		{ // for scoping of the counter.
			auto counter = 0;
			for (auto iter=functions_.begin(); iter!=functions_.end(); iter++, counter++) {
				value(counter) = (*iter)->Eval<T>();
			}
		}

		return value;
	}


	
	// these two lines are explicit instantiations of the template above.  template definitions separate from declarations cause linking problems.  
	// see
	// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
	template Vec<dbl> System::Eval(const Vec<dbl> & variable_values, dbl const& path_variable_value);
	template Vec<mpfr> System::Eval(const Vec<mpfr> & variable_values, mpfr const& path_variable_value);



	template<typename T>
	Mat<T> System::Jacobian(const Vec<T> & variable_values)
	{
		if (variable_values.size()!=NumVariables())
			throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

		if (have_path_variable_)
			throw std::runtime_error("not using a time value for computation of jacobian, but a path variable is defined.");


		if (!is_differentiated_)
		{
			jacobian_.resize(NumFunctions());
			for (int ii = 0; ii < NumFunctions(); ++ii)
			{
				jacobian_[ii] = std::make_shared<bertini::node::Jacobian>(functions_[ii]->Differentiate());
			}
			is_differentiated_ = true;
		}

		SetVariables(variable_values);

		auto vars = Variables(); //TODO: replace this with something that peeks directly into the variables without this copy.

		Mat<T> J(NumFunctions(), NumVariables());
		for (int ii = 0; ii < NumFunctions(); ++ii)
			for (int jj = 0; jj < NumVariables(); ++jj)
				J(ii,jj) = jacobian_[ii]->EvalJ<T>(vars[jj]);

		return J;
	}

	// these two lines are explicit instantiations of the template above.  template definitions separate from declarations cause linking problems.  
	// see
	// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
	template Mat<dbl> System::Jacobian(const Vec<dbl> & variable_values);
	template Mat<mpfr> System::Jacobian(const Vec<mpfr> & variable_values);

	template<typename T>
	Mat<T> System::Jacobian(const Vec<T> & variable_values, const T & path_variable_value)
	{
		if (variable_values.size()!=NumVariables())
			throw std::runtime_error("trying to evaluate jacobian, but number of variables doesn't match.");

		if (!have_path_variable_)
			throw std::runtime_error("trying to use a time value for computation of jacobian, but no path variable defined.");


		if (!is_differentiated_)
		{
			jacobian_.resize(NumFunctions());
			for (int ii = 0; ii < NumFunctions(); ++ii)
			{
				jacobian_[ii] = std::make_shared<bertini::node::Jacobian>(functions_[ii]->Differentiate());
			}
			is_differentiated_ = true;
		}

		SetVariables(variable_values);
		SetPathVariable(path_variable_value);

		auto vars = Variables(); //TODO: replace this with something that peeks directly into the variables without this copy.
		
		Mat<T> J(NumFunctions(), NumVariables());
		for (int ii = 0; ii < NumFunctions(); ++ii)
			for (int jj = 0; jj < NumVariables(); ++jj)
				J(ii,jj) = jacobian_[ii]->EvalJ<T>(vars[jj]);

		return J;

	}

	// these two lines are explicit instantiations of the template above.  template definitions separate from declarations cause linking problems.  
	// see
	// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
	template Mat<dbl> System::Jacobian(const Vec<dbl> & variable_values, const dbl & path_variable_value);
	template Mat<mpfr> System::Jacobian(const Vec<mpfr> & variable_values, const mpfr & path_variable_value);




	template<typename T>
	Vec<T> System::TimeDerivative(const Vec<T> & variable_values, const T & path_variable_value)
	{
		if (!HavePathVariable())
			throw std::runtime_error("computing time derivative of system with no path variable defined");

		if (!is_differentiated_)
		{
			jacobian_.resize(NumFunctions());
			for (int ii = 0; ii < NumFunctions(); ++ii)
			{
				jacobian_[ii] = std::make_shared<bertini::node::Jacobian>(functions_[ii]->Differentiate());
			}
			is_differentiated_ = true;
		}


		SetVariables(variable_values);
		SetPathVariable(path_variable_value);

		Vec<T> ds_dt(NumFunctions());


		for (int ii = 0; ii < NumFunctions(); ++ii)
			ds_dt(ii) = jacobian_[ii]->EvalJ<T>(path_variable_);		


		return ds_dt;
	}


	// these two lines are explicit instantiations of the template above.  template definitions separate from declarations cause linking problems.  
	// see
	// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
	template Vec<dbl> System::TimeDerivative(const Vec<dbl> & variable_values, const dbl & path_variable_value);
	template Vec<mpfr> System::TimeDerivative(const Vec<mpfr> & variable_values, const mpfr & path_variable_value);




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
		for (auto curr_function : functions_)
		{	
			for (auto curr_var_gp : hom_variable_groups_)
			{
				if (!curr_function->IsHomogeneous(curr_var_gp))
					throw std::runtime_error("inhomogeneous function, with homogeneous variable group");

				if (!curr_function->IsPolynomial(curr_var_gp))
					throw std::runtime_error("non-polynomial function, with homogeneous variable group");
			}

			for (auto curr_var_gp : variable_groups_)
			{
				if (!curr_function->IsPolynomial(curr_var_gp))
					throw std::runtime_error("non-polynomial function, with homogeneous variable group");
			}
		}

		bool already_had_homvars = NumHomVariables()!=0;
		
		if (already_had_homvars && NumHomVariables()!=NumVariableGroups())
			throw std::runtime_error("size mismatch on number of homogenizing variables and number of variable groups");


		if (!already_had_homvars)
		{
			homogenizing_variables_.resize(NumVariableGroups());
			// already_had_homvars = false;
		}

		auto group_counter = 0;
		for (auto curr_var_gp = variable_groups_.begin(); curr_var_gp!=variable_groups_.end(); curr_var_gp++)
		{
			std::stringstream converter;
			converter << "HOM_VAR_" << group_counter;

			Var hom_var;
			if (already_had_homvars){
				hom_var = homogenizing_variables_[group_counter];
				if (hom_var->name()==converter.str())
					throw std::runtime_error("duplicate names for homogenizing variables in system during homogenization");
			}
			else
			{
				hom_var = std::make_shared<bertini::node::Variable>(converter.str());
				homogenizing_variables_[group_counter] = hom_var;
			}
		
			for (auto curr_function : functions_)
				curr_function->Homogenize(*curr_var_gp, hom_var);

			group_counter++;
		}

		// need to add patch equations for all variable groups

	}




	bool System::IsHomogeneous() const
	{
		bool have_homvars = NumHomVariables()!=0;
		if (have_homvars && NumHomVariables()!=NumVariableGroups())
			throw std::runtime_error("trying to check homogoneity on a partially-formed system.  mismatch between number of homogenizing variables, and number of variable groups");


		for (auto iter : functions_)
		{
			auto counter = 0;
			for (auto vars : variable_groups_)
			{
				auto tempvars = vars;
				if (have_homvars)
					tempvars.push_front(homogenizing_variables_[counter]);
				counter++;

				if (!iter->IsHomogeneous(tempvars))
					return false;

			}
			for (auto vars : hom_variable_groups_)
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


		for (auto iter : functions_)
		{
			auto counter = 0;
			for (auto vars : variable_groups_)
			{
				auto tempvars = vars;
				if (have_homvars)
					tempvars.push_front(homogenizing_variables_[counter]);

				counter++;

				if (!iter->IsPolynomial(tempvars))
					return false;

			}
			for (auto vars : hom_variable_groups_)
				if (!iter->IsPolynomial(vars))
					return false;

		}
		return true;
	}
	



	//////////////////
	//
	// templated setters
	//
	/////////////////////////////


	template<typename T>
	void System::SetVariables(const Vec<T> & new_values)
	{
		assert(new_values.size()== NumVariables());

		auto vars = Variables();

		auto counter = 0;

		for (auto iter=vars.begin(); iter!=vars.end(); iter++, counter++) {
			(*iter)->set_current_value(new_values(counter));
		}

	}


	// these two lines are explicit instantiations of the template above.  template definitions separate from declarations cause linking problems.  
	// see
	// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
	template void System::SetVariables(const Vec<dbl> & new_values);
	template void System::SetVariables(const Vec<mpfr> & new_values);


	template<typename T>
	void System::SetPathVariable(T new_value)
	{
		if (!have_path_variable_)
			throw std::runtime_error("trying to set the value of the path variable, but one is not defined for this system");

		path_variable_->set_current_value(new_value);
	}



	// these two lines are explicit instantiations of the template above.  template definitions separate from declarations cause linking problems.  
	// see
	// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
	template void System::SetPathVariable(dbl T);
	template void System::SetPathVariable(mpfr T);


	template<typename T>
	void System::SetImplicitParameters(Vec<T> new_values)
	{
		if (new_values.size()!= implicit_parameters_.size())
			throw std::runtime_error("trying to set implicit parameter values, but there is a size mismatch");

		size_t counter = 0;
		for (auto iter=implicit_parameters_.begin(); iter!=implicit_parameters_.end(); iter++, counter++)
			(*iter)->set_current_value(new_values(counter));

	}

	template void System::SetImplicitParameters(Vec<dbl> new_values);
	template void System::SetImplicitParameters(Vec<mpfr> new_values);








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
		time_order_of_variable_groups_.push_back( VariableGroupType::Affine);
	}




	void System::AddHomVariableGroup(VariableGroup const& v)
	{
		hom_variable_groups_.push_back(v);
		is_differentiated_ = false;
		have_ordering_ = false;
		time_order_of_variable_groups_.push_back( VariableGroupType::Homogeneous);
	}





	void System::AddUngroupedVariable(Var const& v)
	{
		ungrouped_variables_.push_back(v);
		is_differentiated_ = false;
		have_ordering_ = false;
		time_order_of_variable_groups_.push_back( VariableGroupType::Ungrouped);
	}




	void System::AddUngroupedVariables(VariableGroup const& v)
	{
		ungrouped_variables_.insert( ungrouped_variables_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
		have_ordering_ = false;
		for (auto iter : v)
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
		Fn F = std::make_shared<node::Function>(N);
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


	VariableGroup System::FIFOVariableOrdering() const
	{
		if (NumHomVariables() != NumVariableGroups())
			throw std::runtime_error("mismatch between number of homogenizing variables, and number of affine variables groups.  unable to form variable vector in FIFO ordering.  you probably need to homogenize the system.");

		VariableGroup ordering;

		unsigned affine_group_counter = 0, ungrouped_variable_counter = 0, hom_group_counter = 0;
		for (auto group_type : time_order_of_variable_groups_)
		{
			switch (group_type){
				case VariableGroupType::Affine:
				{
					ordering.push_back(homogenizing_variables_[affine_group_counter]);
					ordering.insert(ordering.end(), 
					                variable_groups_[affine_group_counter].begin(), 
					                variable_groups_[affine_group_counter].end());
					affine_group_counter++;
				}
				case VariableGroupType::Homogeneous:
				{
					ordering.insert(ordering.end(), 
					                hom_variable_groups_[hom_group_counter].begin(), 
					                hom_variable_groups_[hom_group_counter].end());
					hom_group_counter++;
				}
				case VariableGroupType::Ungrouped:
				{
					ordering.push_back(ungrouped_variables_[ungrouped_variable_counter]);
					ungrouped_variable_counter++;
				}
				default:
				{
					throw std::runtime_error("unacceptable VariableGroupType in FIFOVariableOrdering");
				}
			}
		}
		
		#ifndef BERTINI_DISABLE_ASSERTS
		assert(ordering.size()==NumVariables());
		#endif

		return ordering;	
	}

	VariableGroup System::CanonicalVariableOrdering() const
	{
		if (NumHomVariables() != NumVariableGroups())
			throw std::runtime_error("mismatch between number of homogenizing variables, and number of affine variables groups.  unable to form variable vector in FIFO ordering.  you probably need to homogenize the system.");

		VariableGroup ordering;

		for (auto var_group=variable_groups_.begin(); var_group!=variable_groups_.end(); var_group++)
		{
			ordering.push_back(*(homogenizing_variables_.begin()+ (var_group-variable_groups_.begin())));
			ordering.insert(ordering.end(),var_group->begin(),var_group->end());
		}

		for (auto var_group=hom_variable_groups_.begin(); var_group!=hom_variable_groups_.end(); var_group++)
			ordering.insert(ordering.end(),var_group->begin(),var_group->end());

		ordering.insert(ordering.end(),ungrouped_variables_.begin(),ungrouped_variables_.end());

		#ifndef BERTINI_DISABLE_ASSERTS
		assert(ordering.size()==NumVariables());
		#endif

	    return ordering;
	}

	
	void System::SetFIFOVariableOrdering()
	{
		variable_ordering_ = FIFOVariableOrdering();
		ordering_ = OrderingChoice::FIFO;
		have_ordering_ = true;
	}

	void System::SetCanonicalVariableOrdering()
	{
		variable_ordering_ = CanonicalVariableOrdering();
		ordering_ = OrderingChoice::Canonical;
		have_ordering_ = true;
	}



	VariableGroup System::Variables() const
	{
		if (!have_ordering_ || ordering_==OrderingChoice::Unset)
			throw std::runtime_error("trying to get variable ordering, but is not set.");

		return variable_ordering_;
	}

	//public
	template<typename T>
    Vec<T> System::DehomogenizePoint(Vec<T> const& x) const
    {

    	if (x.size()!=NumVariables())
    		throw std::runtime_error("dehomogenizing point with incorrect number of coordinates");

    	if (!have_ordering_)
			throw std::runtime_error("dehomogenizing point, but variable ordering is not set.");



		switch (ordering_)
		{
			case OrderingChoice::Canonical:
			{
				return DehomogenizePointCanonicalOrdering(x);
			}
			case OrderingChoice::FIFO:
			{
				return DehomogenizePointFIFOOrdering(x);
			}
			default: //ordering_==OrderingChoice::Unset
			{
				throw std::runtime_error("invalid ordering choice");
			}
		}
    }

    template Vec<dbl> System::DehomogenizePoint(Vec<dbl> const& x) const;
	template Vec<mpfr> System::DehomogenizePoint(Vec<mpfr> const& x) const;



    // private
    template<typename T>
    Vec<T> System::DehomogenizePointCanonicalOrdering(Vec<T> const& x) const
    {
    	#ifndef BERTINI_DISABLE_ASSERTS
    	assert(ordering_==OrderingChoice::Canonical);
    	#endif

    	Vec<T> x_dehomogenized(NumNaturalVariables());

    	unsigned hom_index = 0; // index into x, the point we are dehomogenizing
    	unsigned dehom_index = 0; // index into x_dehomogenized, the point we are computing

    	for (auto iter=variable_groups_.begin(); iter!=variable_groups_.end(); iter++)
		{
			auto h = x(hom_index++);
			for (unsigned ii = 0; ii < iter->size(); ++ii)
				x_dehomogenized(dehom_index++) = x(hom_index++) / h;
		}

		for (auto iter=hom_variable_groups_.begin(); iter!=hom_variable_groups_.end(); iter++)
			for (unsigned ii = 0; ii < iter->size(); ++ii)
				x_dehomogenized(dehom_index++) = x(hom_index++);
		
		for (unsigned ii = 0; ii < ungrouped_variables_.size(); ++ii)
			x_dehomogenized(dehom_index++) = x(hom_index++);

		return x_dehomogenized;
    }
    template Vec<dbl> System::DehomogenizePointCanonicalOrdering(Vec<dbl> const& new_values) const;
	template Vec<mpfr> System::DehomogenizePointCanonicalOrdering(Vec<mpfr> const& new_values) const;


    // private
    template<typename T>
    Vec<T> System::DehomogenizePointFIFOOrdering(Vec<T> const& x) const
    {
    	#ifndef BERTINI_DISABLE_ASSERTS
    	assert(ordering_==OrderingChoice::FIFO);
    	#endif

    	Vec<T> x_dehomogenized(NumNaturalVariables());

    	unsigned affine_group_counter = 0;
    	unsigned ungrouped_variable_counter = 0;

    	unsigned hom_index = 0; // index into x, the point we are dehomogenizing
    	unsigned dehom_index = 0; // index into x_dehomogenized, the point we are computing



		for (auto iter : time_order_of_variable_groups_)
		{
			switch (iter){
				case VariableGroupType::Affine:
				{
					// ordering.push_back(homogenizing_variables_[affine_group_counter]);
					// ordering.insert(ordering.end(), iter.begin(), iter.end());
					affine_group_counter++;
				}
				case VariableGroupType::Homogeneous:
				{
					// ordering.insert(ordering.end(), iter.begin(), iter.end());
				}
				case VariableGroupType::Ungrouped:
				{
					// ordering.push_back(ungrouped_variables_[ungrouped_variable_counter]);
					ungrouped_variable_counter++;
				}
				default:
				{
					throw std::runtime_error("unacceptable VariableGroupType in FIFOVariableOrdering");
				}
			}
		}


		return x_dehomogenized;
    }
    template Vec<dbl> System::DehomogenizePointFIFOOrdering(Vec<dbl> const& new_values) const;
	template Vec<mpfr> System::DehomogenizePointFIFOOrdering(Vec<mpfr> const& new_values) const;



	std::vector<int> System::Degrees() const
	{
		std::vector<int> deg;
		for (auto iter : functions_)
		{
			deg.push_back(iter->Degree());
		}
		return deg;
	}


	std::vector<int> System::Degrees(VariableGroup const& vars) const
	{
		std::vector<int> deg;
		for (auto iter : functions_)
		{
			deg.push_back(iter->Degree(vars));
		}
		return deg;
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



	void System::ClearVariables()
	{
		ungrouped_variables_.clear();
		variable_groups_.clear();
		hom_variable_groups_.clear();
		homogenizing_variables_.clear();

		path_variable_.reset();
		have_path_variable_ = false;
	}



	void System::CopyVariableStructure(System const& other)
	{
		this->ClearVariables();

		this->ungrouped_variables_ = other.ungrouped_variables_;
		this->variable_groups_ = other.variable_groups_;
		hom_variable_groups_ = other.hom_variable_groups_;
		homogenizing_variables_ = other.homogenizing_variables_;

		path_variable_ = other.path_variable_;
		have_path_variable_ = other.have_path_variable_;

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
		for (auto iter : s.variable_groups_)
		{
			out << "group " << counter << ": "<< "\n";
			for (auto jter : iter)
				out << *jter << " ";


			out << "\n";
			counter++;
		}

		out << s.NumHomVariables() << " homogenizing variables:\n";
		for (auto iter : s.homogenizing_variables_)
			out << iter << " ";
		out << "\n";


		out << s.NumFunctions() << " functions:\n";
		for (auto iter : s.functions_) 
			out << (iter)->name() << " = " << *iter << "\n";
		out << "\n";


		if (s.NumParameters()) {
			out << s.NumParameters() << " explicit parameters:\n";
			for (auto iter : s.explicit_parameters_)
				out << (iter)->name() << " = " << *iter << "\n";
			out << "\n";
		}


		if (s.NumConstants()) {
			out << s.NumConstants() << " constants:\n";
			for (auto iter : s.constant_subfunctions_)
				out << (iter)->name() << " = " << *iter << "\n";
			out << "\n";
		}

		if (s.path_variable_)
			out << "path variable defined.  named " << s.path_variable_->name() << "\n";

		if (s.is_differentiated_)
		{
			for (auto iter : s.jacobian_) {
				out << (iter)->name() << " = " << *iter << "\n";
			}
			out << "\n";
		}
		else{
			out << "system not differentiated\n";
		}


		return out;
	}







	System System::operator+=(System const& rhs)
	{
		if (this->NumFunctions()!=rhs.NumFunctions())
			throw std::runtime_error("cannot add two Systems with differing numbers of functions");

		if (this->NumVariables()!=rhs.NumVariables())
			throw std::runtime_error("cannot add two Systems with differing numbers of variables");

		for (auto iter=functions_.begin(); iter!=functions_.end(); iter++)
		{
			(*iter)->SetRoot( (*(rhs.functions_.begin()+(iter-functions_.begin())))->entry_node() + (*iter)->entry_node());
		}

		return *this;
	}

	System operator+(System lhs, System const& rhs)
	{
		return lhs+=rhs;
	}


	System System::operator*=(std::shared_ptr<node::Node> const& N)
	{
		for (auto iter=functions_.begin(); iter!=functions_.end(); iter++)
		{
			(*iter)->SetRoot( N * (*iter)->entry_node());
		}
		return *this;
	}


	System operator*(System s, std::shared_ptr<node::Node> const&  N)
	{
		return s*=N;
	}


	System operator*(std::shared_ptr<node::Node> const&  N, System const& s)
	{
		return s*N;
	}


	






}
