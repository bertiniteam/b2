#include "system.hpp"


template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

namespace bertini 
{

	void System::precision(unsigned new_precision)
	{
		for (auto iter : functions_) {
//				iter->precision(new_precision);
		}

		for (auto iter : subfunctions_) {
			//				iter->precision(new_precision);
		}

		for (auto iter : explicit_parameters_) {
			//				iter->precision(new_precision);
		}

		for (auto iter : variables_) {
			//				iter->precision(new_precision);
		}

		for (auto iter :implicit_parameters_) {
			//				iter->precision(new_precision);
		}

		for (auto iter : constant_subfunctions_) {
			//				iter->precision(new_precision);
			//				iter->eval<>()
		}

//			path_variable_->precision(new_precision);

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
		for (auto curr_var_gp : variable_groups_)
		{
			std::stringstream converter;
			converter << "HOM_VAR_" << group_counter;
			Var hom_var = std::make_shared<Variable>(converter.str());

			for (auto curr_function : functions_)
			{
				curr_function->Homogenize(curr_var_gp, hom_var);
				curr_var_gp.push_front(hom_var);
			}

			group_counter++;
		}

		// need to add patch equations for all variable groups

	}




	std::ostream& operator<<(std::ostream& out, const bertini::System & s)
	{
		out << "system:\n\n";
		out << s.NumVariables() << " variables:\n";
		for (auto iter : s.variables_) {
			out << (iter)->name() << "\n";
		}
		out << "\n";

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
