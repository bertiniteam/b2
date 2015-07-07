#include "system.hpp"


template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

namespace bertini 
{



	/////////
	//
	//  getters
	//
	/////////////////


	auto System::NumFunctions() const
	{
		return functions_.size();
	}


	auto System::NumVariables() const
	{
		return variables_.size();
	}


	auto System::NumVariableGroups() const
	{
		return variable_groups_.size();
	}

	auto System::NumHomVariableGroups() const
	{
		return hom_variable_groups_.size();
	}


	auto System::NumConstants() const
	{
		return constant_subfunctions_.size();
	}

	auto System::NumParameters() const
	{
		return explicit_parameters_.size();
	}


	/**
	 Get the number of implicit parameters in this system
	 */
	auto System::NumImplicitParameters() const
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

		for (auto iter : variables_) {
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


		for (auto iter : ungrouped_variables_)
			iter->precision(new_precision);

		for (auto iter : variable_groups_)
			for (auto jter : iter)
				jter->precision(new_precision);

		for (auto iter : hom_variable_groups_)
			for (auto jter : iter)
				jter->precision(new_precision);

		precision_ = new_precision;
	}



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
				jacobian_[ii] = std::make_shared<bertini::Jacobian>(functions_[ii]->Differentiate());
			}
			is_differentiated_ = true;
		}

		SetVariables(variable_values);

		Mat<T> J(NumFunctions(), NumVariables());
		for (int ii = 0; ii < NumFunctions(); ++ii)
		{
			for (int jj = 0; jj < NumVariables(); ++jj)
			{
				J(ii,jj) = jacobian_[ii]->EvalJ<T>(variables_[jj]);
			}
		}

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
				jacobian_[ii] = std::make_shared<bertini::Jacobian>(functions_[ii]->Differentiate());
			}
			is_differentiated_ = true;
		}

		SetVariables(variable_values);
		SetPathVariable(path_variable_value);

		Mat<T> J(NumFunctions(), NumVariables());
		for (int ii = 0; ii < NumFunctions(); ++ii)
		{
			for (int jj = 0; jj < NumVariables(); ++jj)
			{
				J(ii,jj) = jacobian_[ii]->EvalJ<T>(variables_[jj]);
			}
		}

		return J;

	}

	// these two lines are explicit instantiations of the template above.  template definitions separate from declarations cause linking problems.  
	// see
	// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
	template Mat<dbl> System::Jacobian(const Vec<dbl> & variable_values, const dbl & path_variable_value);
	template Mat<mpfr> System::Jacobian(const Vec<mpfr> & variable_values, const mpfr & path_variable_value);

	void System::Homogenize()
	{

		
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
		


		auto group_counter = 0;
		for (auto curr_var_gp = variable_groups_.begin(); curr_var_gp!=variable_groups_.end(); curr_var_gp++)
		{
			std::stringstream converter;
			converter << "HOM_VAR_" << group_counter;
			Var hom_var = std::make_shared<Variable>(converter.str());

			for (auto curr_function : functions_)
			{
				curr_function->Homogenize(*curr_var_gp, hom_var);
				
			}
			curr_var_gp->push_front(hom_var);
			variables_.push_back(hom_var);
			group_counter++;
		}

		// need to add patch equations for all variable groups

	}




	bool System::IsHomogeneous() const
	{
		for (auto iter : functions_)
		{
			if (!iter->IsHomogeneous(variables_))
			{
				return false;
			}
			for (auto vars : variable_groups_)
			{
				if (!iter->IsHomogeneous(vars))
				{
					return false;
				}
			}
			for (auto vars : hom_variable_groups_)
			{
				if (!iter->IsHomogeneous(vars))
				{
					return false;
				}
			}

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
		assert(new_values.size()== variables_.size());

		{ // for scoping of the counter.
			auto counter = 0;
			for (auto iter=variables_.begin(); iter!=variables_.end(); iter++, counter++) {
				(*iter)->set_current_value(new_values(counter));
			}
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
		path_variable_->set_current_value(new_value);
		have_path_variable_ = true;
	}



	// these two lines are explicit instantiations of the template above.  template definitions separate from declarations cause linking problems.  
	// see
	// https://isocpp.org/wiki/faq/templates#templates-defn-vs-decl
	template void System::SetPathVariable(dbl T);
	template void System::SetPathVariable(mpfr T);


	template<typename T>
	void System::SetImplicitParameters(Vec<T> new_values)
	{
		assert(new_values.size()== implicit_parameters_.size());

		{
			auto counter = 0;
			for (auto iter=implicit_parameters_.begin(); iter!=implicit_parameters_.end(); iter++, counter++) {
				(*iter)->set_current_value(new_values(counter));
			}
		}
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
		variables_.insert( variables_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
	}




	void System::AddHomVariableGroup(VariableGroup const& v)
	{
		hom_variable_groups_.push_back(v);
		variables_.insert( variables_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
	}





	void System::AddUngroupedVariable(Var const& v)
	{
		ungrouped_variables_.push_back(v);
		variables_.push_back(v);
		is_differentiated_ = false;
	}




	void System::AddUngroupedVariables(VariableGroup const& v)
	{
		ungrouped_variables_.insert( ungrouped_variables_.end(), v.begin(), v.end() );
		variables_.insert( variables_.end(), v.begin(), v.end() );
		is_differentiated_ = false;
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
		Fn F = std::make_shared<Function>(N);
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








	std::vector<int> System::Degrees() const
	{
		std::vector<int> deg;
		for (auto iter : functions_)
		{
			deg.push_back(iter->Degree());
		}
		return deg;
	}



	void System::ReorderFunctionsByDegreeDecreasing()
	{
		auto degs = Degrees();

		// now we sort a vector of the indexing numbers by the degrees contained in degs.
		std::vector<size_t> indices(degs.size());
		//http://en.cppreference.com/w/cpp/algorithm/iota
		//http://www.cplusplus.com/doc/tutorial/typecasting/
		std::iota(begin(indices), end(indices), static_cast<size_t>(0));
		std::sort( begin(indices), end(indices), [&](size_t a, size_t b) { return degs[a] > degs[b]; } );
		


		// finally, we re-order the functions based on the indices we just computed
		std::vector<std::shared_ptr<Function> > re_ordered_functions(degs.size());
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
		auto degs = Degrees();

		// now we sort a vector of the indexing numbers by the degrees contained in degs.
		std::vector<size_t> indices(degs.size());
		//http://en.cppreference.com/w/cpp/algorithm/iota
		//http://www.cplusplus.com/doc/tutorial/typecasting/
		std::iota(begin(indices), end(indices), static_cast<size_t>(0));
		std::sort( begin(indices), end(indices), [&](size_t a, size_t b) { return degs[a] < degs[b]; } );
		


		// finally, we re-order the functions based on the indices we just computed
		std::vector<std::shared_ptr<Function> > re_ordered_functions(degs.size());
		size_t ind = 0;
		for (auto iter : indices)
		{
			re_ordered_functions[ind] = functions_[iter];
			ind++;
		}

		swap(functions_, re_ordered_functions);
	}








	
	//////////////////
	//
	//  output operators
	//
	////////////////////

	std::ostream& operator<<(std::ostream& out, const bertini::System & s)
	{
		out << "system:\n\n";
		out << s.NumVariables() << " variables:\n";
		for (auto iter : s.variables_) {
			out << (iter)->name() << "\n";
		}
		out << "\n";


		out << s.NumVariableGroups() << " variable groups, containing these variables:\n";
		auto counter = 0;
		for (auto iter : s.variable_groups_)
		{
			out << "group" << counter << " - "<< "\n";
			for (auto jter : iter)
			{
				out << *jter << " ";
			}
			out << "\n";
			counter++;
		}

		out << s.NumFunctions() << " functions:\n";
		for (auto iter : s.functions_) {
			out << (iter)->name() << "\n";
			out << *iter << "\n";
		}
		out << "\n";


		if (s.NumParameters()) {
			out << s.NumParameters() << " explicit parameters:\n";
			for (auto iter : s.explicit_parameters_) {
				out << (iter)->name() << "\n";
				out << *iter << "\n";
			}
			out << "\n";
		}


		if (s.NumConstants()) {
			out << s.NumConstants() << " constants:\n";
			for (auto iter : s.constant_subfunctions_) {
				out << (iter)->name() << "\n";
				out << *iter << "\n";
			}
			out << "\n";
		}

		if (s.path_variable_) {
			out << "path variable defined.  named " << s.path_variable_->name() << "\n";
		}

		if (s.is_differentiated_)
		{
			for (auto iter : s.jacobian_) {
				out << (iter)->name() << "\n";
				out << *iter << "\n";
			}
			out << "\n";
		}
		else{
			out << "system not differentiated\n";
		}


		return out;
	}
}
