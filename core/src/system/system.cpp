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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire


#include "bertini2/system/system.hpp"

template<typename NumT> using Vec = bertini::Vec<NumT>;
template<typename NumT> using Mat = bertini::Mat<NumT>;
using Nd = std::shared_ptr<bertini::node::Node>;

BOOST_CLASS_EXPORT(bertini::System)



namespace bertini 
{

	using namespace bertini::node;
	
	EvalMethod DefaultEvalMethod()
	{
		return EvalMethod::SLP;
	}

	DerivMethod DefaultDerivMethod()
	{
		return DerivMethod::Derivatives;
	}

	bool DefaultAutoSimplify()
	{
		// differentiation now emits already-simplified trees, so the post-hoc
		// Simplify pass is redundant -- and it mutates IN PLACE, including
		// subtrees the derivative shares with the user's original functions.
		// holding f must never observe f changing because the system was
		// differentiated.  Simplify()/AutoSimplify(true) remain explicit opt-ins.
		return false;
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

		swap(a.eval_method_,b.eval_method_);

		swap(a.precision_,b.precision_);
		swap(a.is_patched_,b.is_patched_);
		swap(a.patch_,b.patch_);
	}

	// construct from a list of functions, auto-discovering the variables
	System::System(std::vector<Fn> const& functions) : System()
	{
		AddFunctions(functions);
		AddVariableGroup( node::GatherVariables(functions) );
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

		eval_method_ = other.eval_method_;

		time_order_of_variable_groups_ = other.time_order_of_variable_groups_;

		current_variable_values_ = other.current_variable_values_;

		variable_ordering_ = other.variable_ordering_;
		have_ordering_ =  other.have_ordering_;

		precision_ = other.precision_;


		constant_subfunctions_ = other.constant_subfunctions_;
		subfunctions_ = other .subfunctions_;
		functions_ = other .functions_;
		explicit_parameters_  = other .explicit_parameters_;

		// // now to do the members which are not simply copied
		// constant_subfunctions_.resize(other.constant_subfunctions_.size());
		// for (unsigned ii = 0; ii < constant_subfunctions_.size(); ++ii)
		// 	constant_subfunctions_[ii] = Function::Make(other.constant_subfunctions_[ii]->EntryNode());

		// subfunctions_.resize(other.subfunctions_.size());
		// for (unsigned ii = 0; ii < subfunctions_.size(); ++ii)
		// 	subfunctions_[ii] = Function::Make(other.subfunctions_[ii]->EntryNode());

		// functions_.resize(other.functions_.size());
		// for (unsigned ii = 0; ii < functions_.size(); ++ii)
		// 	functions_[ii] = Function::Make(other.functions_[ii]->EntryNode());

		// explicit_parameters_.resize(other.explicit_parameters_.size());
		// for (unsigned ii = 0; ii < explicit_parameters_.size(); ++ii)
		// 	explicit_parameters_[ii] = Function::Make(other.explicit_parameters_[ii]->EntryNode());
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


	size_t System::NumNaturalFunctions() const
	{
		if (!blocks_.empty())
		{
			size_t n = 0;
			for (auto const& blk : blocks_)
				n += std::visit([](auto const& b){ return b.NumFunctions(); }, blk);
			return n;
		}
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
		return NumNaturalFunctions() + NumPatches();
	}


	void System::precision(unsigned new_precision) const
	{
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

		if (!blocks_.empty())
		{
			for (auto const& blk : blocks_)
				std::visit([&](auto const& b){ b.Precision(new_precision); }, blk);
		}
		else
		switch (eval_method_)
		{
			case EvalMethod::FunctionTree:{

				if (is_differentiated_)
				{
					switch (deriv_method_){
						case DerivMethod::JacobianNode:{
							for (const auto& iter : jacobian_)
								iter->precision(new_precision);
							break;
						}
						case DerivMethod::Derivatives:{
							for (const auto& iter : space_derivatives_)
								iter->precision(new_precision);
							for (const auto& iter : time_derivatives_)
								iter->precision(new_precision);
							break;
						}
					}
				}
				break;
			}
			case EvalMethod::SLP:
			{
				// the SLP exists (and is used for plain Eval) regardless of whether
				// the system has been differentiated, so its precision must be kept
				// in sync unconditionally.  previously this was gated behind
				// is_differentiated_, leaving a never-differentiated system's SLP at
				// its compile-time precision forever.
				this->slp_.precision(new_precision);
				break;
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
		Precision(std::get<Vec<mpfr_complex> >(current_variable_values_),new_precision);

		if (IsPatched())
			patch_.Precision(new_precision);

		precision_ = new_precision;
	}


	void System::Differentiate() const
	{
		if (!blocks_.empty()) { is_differentiated_ = true; return; }

		switch (deriv_method_){
			case DerivMethod::JacobianNode:
			{
				DifferentiateUsingJacobianNode();
				break;
			}
			case DerivMethod::Derivatives:
			{
				DifferentiateUsingDerivatives();
				break;
			}
		}


		if (auto_simplify_)
			this->SimplifyDerivatives();


		switch (eval_method_)
		{
			case EvalMethod::FunctionTree:{
				break;
			}
			case EvalMethod::SLP:
			{	
				SLPCompiler compiler;
				this->slp_ = compiler.Compile(*this);
				break;
			}
		}
		


	}

	void System::DifferentiateUsingJacobianNode() const
	{
		auto num_functions = NumNaturalFunctions();
		jacobian_.resize(num_functions);
		for (size_t ii = 0; ii < num_functions; ++ii)
			jacobian_[ii] = Jacobian::Make(functions_[ii]->Differentiate());

		is_differentiated_ = true;
	}

	void System::DifferentiateUsingDerivatives() const
	{
		const auto& vars = this->Variables();
		const auto num_vars = NumVariables();
		const auto num_functions = NumNaturalFunctions();

		space_derivatives_.resize(num_functions*num_vars);
		// again, computing these in column major, so staying with one variable at a time.
		for (size_t jj = 0; jj < num_vars; ++jj)
			for (size_t ii = 0; ii < num_functions; ++ii)
				space_derivatives_[ii+jj*num_functions] = Function::Make(functions_[ii]->Differentiate(vars[jj]));

		if (HavePathVariable())
		{
			const auto& t = path_variable_;
			time_derivatives_.resize(num_functions);
				for (size_t ii = 0; ii < num_functions; ++ii)
					time_derivatives_[ii] = Function::Make(functions_[ii]->Differentiate(t));
		}

		is_differentiated_ = true;
	}

	std::vector< Nd > System::GetSpaceDerivatives() const
	{
		if ( (deriv_method_==DerivMethod::JacobianNode) || (!is_differentiated_) )
			DifferentiateUsingDerivatives();

		return space_derivatives_;
	}

	std::vector< Nd > System::GetTimeDerivatives() const
	{
		if ( (deriv_method_==DerivMethod::JacobianNode) || (!is_differentiated_) )
			DifferentiateUsingDerivatives();
		
		return time_derivatives_;
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

		// idempotency: homogenizing an already-homogenized system must be a no-op.
		// without this, a second call re-homogenizes each function with respect to
		// the group INCLUDING its homogenizing variable, inflating degrees (observed
		// 2026-06-06: degrees (2,1) -> (3,2), turning a 2-path TD into a 6-path one).
		if (already_had_homvars && IsHomogeneous())
			return;

		if (!already_had_homvars)
		{
			homogenizing_variables_.resize(NumVariableGroups());
		}




		auto PushFront = [&](auto & container, auto item){
			container.push_back(item);
			std::rotate(container.rbegin(), container.rbegin() + 1, container.rend());
		};




		auto group_counter = 0;
		for (auto curr_var_gp = variable_groups_.begin(); curr_var_gp!=variable_groups_.end(); curr_var_gp++)
		{
			std::stringstream converter;
			converter << "HOM_VAR_" << group_counter;

			if (already_had_homvars){
				Var hom_var = homogenizing_variables_[group_counter];
				VariableGroup temp_group = *curr_var_gp;

				PushFront(temp_group, hom_var);

				// temp_group.push_front(hom_var);
				for (const auto& curr_function : functions_)
					curr_function->Homogenize(temp_group, hom_var);
			}
			else
			{
				Var hom_var = Variable::Make(converter.str());
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
		auto PushFront = [&](auto & container, auto item){
			container.push_back(item);
			std::rotate(container.rbegin(), container.rbegin() + 1, container.rend());
		};


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
					PushFront(tempvars, homogenizing_variables_[counter]);
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

		auto PushFront = [&](auto & container, auto item){
			container.push_back(item);
			std::rotate(container.rbegin(), container.rbegin() + 1, container.rend());
		};



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
					PushFront(tempvars,homogenizing_variables_[counter]);

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




	void System::SetVariableGroups(std::vector<VariableGroup> const& groups)
	{
		// clear the existing variable structure, but preserve the path variable.
		ungrouped_variables_.clear();
		variable_groups_.clear();
		hom_variable_groups_.clear();
		homogenizing_variables_.clear();
		time_order_of_variable_groups_.clear();

		// install the supplied groups as affine variable groups.  AddVariableGroup
		// takes care of the FIFO time-ordering entries and resets the relevant flags.
		for (auto const& g : groups)
			AddVariableGroup(g);

		is_differentiated_ = false;
		have_ordering_ = false;
		is_patched_ = false;
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
		time_order_of_variable_groups_.insert(time_order_of_variable_groups_.end(), v.size(), VariableGroupType::Ungrouped);
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



	void System::AddFunction(Nd const& N, std::string const& name)
	{
		functions_.push_back(Function::Make(N, name));
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
		if (constructed_ordering.size()!=NumVariables()){
			std::stringstream err_msg;
			err_msg << "resulting constructed ordering has differing size (" << constructed_ordering.size() << ") from the number of variables in the problem (" << NumVariables() << ").";
			throw std::runtime_error(err_msg.str());
		}

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

		unsigned hom_group_counter(0), affine_group_counter(0);
		for (auto curr_grouptype : time_order_of_variable_groups_)
		{
			if (curr_grouptype==VariableGroupType::Homogeneous)
				s.push_back(static_cast<unsigned>(hom_variable_groups_[hom_group_counter++].size()));
			else if (curr_grouptype==VariableGroupType::Affine)
				s.push_back(static_cast<unsigned>(variable_groups_[affine_group_counter++].size() + (have_homvars ? 1 : 0)));
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
		
		using RealT = typename Eigen::NumTraits<NumT>::Real;
		using ComplexT = NumT;

		RealT bound(0);

		for (unsigned ii=0; ii < num_evaluations; ii++)
		{	
			Vec<ComplexT> randy = RandomOfUnits<ComplexT>(static_cast<unsigned>(NumVariables()));
			Vec<ComplexT> f_vals;
			if (HavePathVariable())
				f_vals = Eval(randy, RandomUnit<ComplexT>());
			else
				f_vals = Eval(randy);
			
			Mat<ComplexT> dh_dx = Jacobian<ComplexT>();
			
			bound = max(f_vals.array().abs().maxCoeff(),
						 dh_dx.array().abs().maxCoeff(), bound);
		}
		return bound;
	}

	template double System::CoefficientBound<dbl>(unsigned) const;
	template mpfr_float System::CoefficientBound<mpfr_complex>(unsigned) const;

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
		is_differentiated_ = false;
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
		is_differentiated_ = false;
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

		is_differentiated_ = false;
		have_ordering_ = false;
	}



	bool System::RemoveVariable(Var const& v)
	{
		// remove the n-th time-ordering entry of the given group type, keeping the
		// FIFO ordering consistent with the variable-group containers.
		auto remove_nth_time_order = [this](VariableGroupType t, size_t n)
		{
			size_t count = 0;
			for (auto it = time_order_of_variable_groups_.begin(); it != time_order_of_variable_groups_.end(); ++it)
			{
				if (*it == t)
				{
					if (count == n)
					{
						time_order_of_variable_groups_.erase(it);
						return;
					}
					++count;
				}
			}
		};

		auto did_remove = [this]()
		{
			is_differentiated_ = false;
			have_ordering_ = false;
			is_patched_ = false;
		};

		// search the affine variable groups
		for (size_t gi = 0; gi < variable_groups_.size(); ++gi)
		{
			auto& group = variable_groups_[gi];
			auto it = std::find(group.begin(), group.end(), v);
			if (it != group.end())
			{
				group.erase(it);
				if (group.empty())
				{
					variable_groups_.erase(variable_groups_.begin() + gi);
					remove_nth_time_order(VariableGroupType::Affine, gi);
				}
				did_remove();
				return true;
			}
		}

		// search the homogeneous / projective variable groups
		for (size_t gi = 0; gi < hom_variable_groups_.size(); ++gi)
		{
			auto& group = hom_variable_groups_[gi];
			auto it = std::find(group.begin(), group.end(), v);
			if (it != group.end())
			{
				group.erase(it);
				if (group.empty())
				{
					hom_variable_groups_.erase(hom_variable_groups_.begin() + gi);
					remove_nth_time_order(VariableGroupType::Homogeneous, gi);
				}
				did_remove();
				return true;
			}
		}

		// search the ungrouped variables (each is its own ungrouped time-order entry)
		{
			auto it = std::find(ungrouped_variables_.begin(), ungrouped_variables_.end(), v);
			if (it != ungrouped_variables_.end())
			{
				size_t idx = static_cast<size_t>(std::distance(ungrouped_variables_.begin(), it));
				ungrouped_variables_.erase(it);
				remove_nth_time_order(VariableGroupType::Ungrouped, idx);
				did_remove();
				return true;
			}
		}

		return false;
	}



	void System::SimplifyFunctions()
	{
		using bertini::Simplify;
		for (auto& iter : this->functions_)
			Simplify(iter);

		is_differentiated_ = false;
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



		switch (deriv_method_){
			case DerivMethod::JacobianNode:{
				for (auto& iter : this->jacobian_)
					Simplify(iter);
				break;
			}
			case DerivMethod::Derivatives:{
				for (auto& iter : this->space_derivatives_)
					Simplify(iter);
				for (auto& iter : this->time_derivatives_)
					Simplify(iter);
				break;
			}
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
		out << "\n";

		out << s.NumHomVariables() << " homogenizing variables:\n";
		for (const auto& iter : s.homogenizing_variables_)
			out << (*iter) << " ";
		out << "\n\n";


		out << s.ungrouped_variables_.size() << " ungrouped variables:\n";
		for (const auto& v :s.ungrouped_variables_)
			out << (*v) << " ";
		out << "\n\n";


		out << s.NumNaturalFunctions() << " functions:\n";
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

				switch (s.deriv_method_){
					case DerivMethod::JacobianNode:{
						out << "using the JacobianNode method of differentiation:" << std::endl;

						for (const auto& iter : s.jacobian_)
							out << (iter)->name() << " = " << *iter << "\n";
						break;
					}

					case DerivMethod::Derivatives:{
						out << "using the Derivatives method of differentiation:" << std::endl;

						for (size_t jj = 0; jj < s.NumVariables(); ++jj)
							for (size_t ii = 0; ii < s.NumNaturalFunctions(); ++ii)
							{
								const auto& d = s.space_derivatives_[ii+jj*s.NumNaturalFunctions()];
								out << "jac_space_der(" << ii << "," << jj << ") = " << d << "\n";
							}

						if (s.HavePathVariable())
							for (size_t ii = 0; ii < s.NumNaturalFunctions(); ++ii)
							{
								const auto& d = s.time_derivatives_[ii];
								out << "jac_time_der(" << ii << ") = " << d << "\n";
							}
						break;
					}
				} // switch on deriv method



				if (s.eval_method_ == EvalMethod::SLP)
				{
					out << "since using SLP for evaluation, here's the SLP:" << std::endl;
					out << s.slp_;				
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

		out << "\ncurrent variable values:\n";
		out << std::get< Vec<dbl> > (s.current_variable_values_) << "\n";
		out << std::get< Vec<mpfr_complex> > (s.current_variable_values_) << "\n";

		return out;
	}









	/////////////////
	//
	// Arithemetic operators
	//
	///////////////////


	System& System::operator+=(System const& rhs)
	{
		if (this->NumTotalFunctions()!=rhs.NumTotalFunctions())
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

		// make NEW Function wrappers rather than calling SetRoot on the existing
		// ones: the existing Function nodes are shared_ptrs, SHARED with whatever
		// system this one was (shallowly) copied from.  mutating them in place
		// rewrites that system's functions too — e.g. forming the homotopy
		// (1-t)*target + gamma*t*start used to corrupt both the target and the
		// start system (observed 2026-06-06).
		for (auto iter=functions_.begin(); iter!=functions_.end(); iter++)
			*iter = node::Function::Make(
				(*(rhs.functions_.begin()+(iter-functions_.begin())))->EntryNode() + (*iter)->EntryNode(),
				(*iter)->name());

		is_differentiated_ = false;
		return *this;
	}

	const System operator+(System lhs, System const& rhs)
	{
		return lhs+=rhs;
	}


	System& System::operator*=(std::shared_ptr<node::Node> const& N)
	{
		// new wrappers, not SetRoot — see comment in operator+= above.
		for (auto iter=functions_.begin(); iter!=functions_.end(); iter++)
		{
			*iter = node::Function::Make( N * (*iter)->EntryNode(), (*iter)->name());
		}
		is_differentiated_ = false;
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

		for (unsigned ii(0); ii<sys2.NumNaturalFunctions(); ++ii)
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

		// Rebuild evaluation machinery from the deserialized expression tree rather
		// than trusting the archived copy: the serialized SLP does not survive the
		// round trip faithfully (its time-derivative outputs read stale memory,
		// observed 2026-06-06; root cause in SLP serialization not yet identified).
		// Differentiate() re-derives the derivative trees and recompiles the SLP
		// from the clone's own (verified-exact) tree.
		if (sys_clone.GetEvalMethod() == EvalMethod::SLP)
			sys_clone.Differentiate();

		// Normalize precision across all parts of the clone.  The source system can
		// carry internally-inconsistent precision state (e.g. precision_ says 30 but
		// the SLP is still at its compile-time precision); precision() propagates to
		// every node, derivative, and the SLP.
		sys_clone.precision(sys_clone.precision());

		return sys_clone;
	}


	void Simplify(System & sys)
	{
		sys.Simplify();
	}


	void System::ResetFunctions() const
	{
		// TODO: it has the unfortunate side effect of resetting constant functions, too.
		switch (eval_method_){
		case EvalMethod::FunctionTree:
			for (const auto& iter : functions_)
				iter->Reset();
			break;
		case EvalMethod::SLP:
			break;
		}
	}

	void System::ResetJacobian() const
	{
		switch (eval_method_)
		{
			case EvalMethod::FunctionTree:{
				switch (deriv_method_){
					case DerivMethod::JacobianNode:
						for (const auto& iter : jacobian_)
							iter->Reset();
						break;
					case DerivMethod::Derivatives:
						for (const auto& iter : space_derivatives_)
							iter->Reset();
						break;
				}
				break;
			}
			case EvalMethod::SLP:
				break;
		}
	}

	void System::ResetTimeDerivatives() const
	{
		switch (eval_method_)
		{
			case EvalMethod::FunctionTree:{
				switch (deriv_method_){
					case DerivMethod::JacobianNode:
						for (const auto& iter : jacobian_)
							iter->Reset();
						break;
					case DerivMethod::Derivatives:
						for (const auto& iter : time_derivatives_)
							iter->Reset();
						break;
				}
				break;
			}
			case EvalMethod::SLP:
				break;
		}
	}

	void System::Reset() const
	{
		ResetFunctions();
		ResetJacobian();
		ResetTimeDerivatives();
	}

	const System::Var& System::GetPathVariable() const
	{
		if (this->HavePathVariable())
			return this->path_variable_;
		throw std::runtime_error("trying to get path variable for a system which doesn't have a path variable defined");
	}


	// Explicit instantiation definitions — paired with extern template declarations in system.hpp.

	template void System::EvalInPlace<dbl>(Vec<dbl>&) const;
	template void System::EvalInPlace<mpfr_complex>(Vec<mpfr_complex>&) const;

	template Vec<dbl> System::Eval<dbl>() const;
	template Vec<mpfr_complex> System::Eval<mpfr_complex>() const;

	template void System::JacobianInPlace<dbl>(Mat<dbl>&) const;
	template void System::JacobianInPlace<mpfr_complex>(Mat<mpfr_complex>&) const;

	template Mat<dbl> System::Jacobian<dbl>() const;
	template Mat<mpfr_complex> System::Jacobian<mpfr_complex>() const;

	template Mat<dbl> System::Jacobian<dbl>(const Vec<dbl>&) const;
	template Mat<mpfr_complex> System::Jacobian<mpfr_complex>(const Vec<mpfr_complex>&) const;

	template void System::JacobianInPlace<dbl>(Mat<dbl>&, const Vec<dbl>&) const;
	template void System::JacobianInPlace<mpfr_complex>(Mat<mpfr_complex>&, const Vec<mpfr_complex>&) const;

	template void System::TimeDerivativeInPlace<dbl>(Vec<dbl>&) const;
	template void System::TimeDerivativeInPlace<mpfr_complex>(Vec<mpfr_complex>&) const;

	template Vec<dbl> System::TimeDerivative<dbl>() const;
	template Vec<mpfr_complex> System::TimeDerivative<mpfr_complex>() const;

	template void System::SetVariables<dbl>(const Vec<dbl>&) const;
	template void System::SetVariables<mpfr_complex>(const Vec<mpfr_complex>&) const;

	template void System::SetPathVariable<dbl>(dbl const&) const;
	template void System::SetPathVariable<mpfr_complex>(mpfr_complex const&) const;

	template void System::SetAndReset<dbl>(Vec<dbl> const&, dbl const&) const;
	template void System::SetAndReset<mpfr_complex>(Vec<mpfr_complex> const&, mpfr_complex const&) const;

	template void System::SetAndReset<dbl>(Vec<dbl> const&) const;
	template void System::SetAndReset<mpfr_complex>(Vec<mpfr_complex> const&) const;

}
