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

		// the polynomial path (functions / derivatives / SLP / eval+deriv method) lives in blocks_
		swap(a.blocks_,b.blocks_);
		swap(a.is_differentiated_,b.is_differentiated_);

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

		is_differentiated_ = other.is_differentiated_;

		time_order_of_variable_groups_ = other.time_order_of_variable_groups_;

		current_variable_values_ = other.current_variable_values_;

		variable_ordering_ = other.variable_ordering_;
		have_ordering_ =  other.have_ordering_;

		precision_ = other.precision_;


		explicit_parameters_  = other .explicit_parameters_;

		// the polynomial path (functions / subfunctions / derivatives / SLP) lives in blocks_
		blocks_ = other.blocks_;

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
		size_t n = 0;
		for (auto const& blk : blocks_)
			n += std::visit([](auto const& b){ return b.NumFunctions(); }, blk);
		return n;
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
		return PolyBlockPtr() ? PolyBlockPtr()->NumConstants() : 0;
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
		for (const auto& iter : explicit_parameters_) {
			iter->precision(new_precision);
		}


		for (const auto& iter :implicit_parameters_) {
			iter->precision(new_precision);
		}

		// each block precisions its own functions / subfunctions / derivatives / SLP
		for (auto const& blk : blocks_)
			std::visit([&](auto const& b){ b.Precision(new_precision); }, blk);

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
		// push the System's variable ordering / path variable / auto-simplify into the
		// polynomial block, then let each block differentiate itself (the polynomial block
		// builds its derivatives + compiles its SLP; structured blocks are analytic no-ops).
		SyncPolyBlock();
		for (auto const& blk : blocks_)
			std::visit([](auto const& b){ b.Differentiate(); }, blk);
		is_differentiated_ = true;
	}

	std::vector< Nd > System::GetSpaceDerivatives() const
	{
		SyncPolyBlock();
		if (auto* p = PolyBlockPtr())
			return p->GetSpaceDerivatives();
		return {};
	}

	std::vector< Nd > System::GetTimeDerivatives() const
	{
		SyncPolyBlock();
		if (auto* p = PolyBlockPtr())
			return p->GetTimeDerivatives();
		return {};
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
		for (const auto& curr_var_gp : hom_variable_groups_)
			for (auto const& b : blocks_)
				if (!std::visit([&](auto const& blk){ return blk.IsHomogeneous(curr_var_gp); }, b))
					throw std::runtime_error("inhomogeneous function, with homogeneous variable group");

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

				// every block homogenizes itself w.r.t. this group: the polynomial block walks
				// its trees, structured blocks fold the constant onto the homogenizing variable.
				for (auto& b : blocks_)
					std::visit([&](auto& blk){ blk.Homogenize(temp_group, hom_var); }, b);
			}
			else
			{
				Var hom_var = Variable::Make(converter.str());
				homogenizing_variables_[group_counter] = hom_var;
				for (auto& b : blocks_)
					std::visit([&](auto& blk){ blk.Homogenize(*curr_var_gp, hom_var); }, b);
			}

			group_counter++;
		}

		InvalidateDifferentiation();
		have_ordering_ = false;

		#ifndef BERTINI_DISABLE_ASSERTS
		assert(homogenizing_variables_.size() == variable_groups_.size());
		#endif
	}


	void System::Homogenize(VariableGroup const& provided_hom_vars)
	{
		// Like Homogenize(), but adopt the supplied homogenizing variables (one per affine variable
		// group, in group order) instead of minting fresh ones.  Mirrors Homogenize()'s
		// fresh-system branch exactly; the only difference is where the hom var comes from.
		for (const auto& curr_var_gp : hom_variable_groups_)
			for (auto const& b : blocks_)
				if (!std::visit([&](auto const& blk){ return blk.IsHomogeneous(curr_var_gp); }, b))
					throw std::runtime_error("inhomogeneous function, with homogeneous variable group");

		if (!IsPolynomial())
			throw std::runtime_error("trying to homogenize a non-polynomial system.");

		if (NumHomVariables()!=0)
			throw std::runtime_error("Homogenize(provided homogenizing variables): system is already homogenized.");

		if (provided_hom_vars.size()!=NumVariableGroups())
			throw std::runtime_error("Homogenize(provided homogenizing variables): need exactly one homogenizing variable per affine variable group.");

		homogenizing_variables_.resize(NumVariableGroups());

		auto group_counter = 0;
		for (auto curr_var_gp = variable_groups_.begin(); curr_var_gp!=variable_groups_.end(); curr_var_gp++)
		{
			Var hom_var = provided_hom_vars[group_counter];
			homogenizing_variables_[group_counter] = hom_var;
			for (auto& b : blocks_)
				std::visit([&](auto& blk){ blk.Homogenize(*curr_var_gp, hom_var); }, b);
			group_counter++;
		}

		InvalidateDifferentiation();
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

		auto all_blocks_homogeneous = [&](VariableGroup const& tempvars) -> bool {
			for (auto const& b : blocks_)
				if (!std::visit([&](auto const& blk){ return blk.IsHomogeneous(tempvars); }, b))
					return false;
			return true;
		};

		auto counter = 0;
		for (const auto& vars : variable_groups_)
		{
			auto tempvars = vars;
			if (have_homvars)
				PushFront(tempvars, homogenizing_variables_[counter]);
			counter++;
			if (!all_blocks_homogeneous(tempvars))
				return false;
		}
		for (const auto& vars : hom_variable_groups_)
			if (!all_blocks_homogeneous(vars))
				return false;
		if (NumUngroupedVariables()>0)
			if (!all_blocks_homogeneous(ungrouped_variables_))
				return false;
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

		auto all_blocks_polynomial = [&](VariableGroup const& tempvars) -> bool {
			for (auto const& b : blocks_)
				if (!std::visit([&](auto const& blk){ return blk.IsPolynomial(tempvars); }, b))
					return false;
			return true;
		};

		auto counter = 0;
		for (const auto& vars : variable_groups_)
		{
			auto tempvars = vars;
			if (have_homvars)
				PushFront(tempvars,homogenizing_variables_[counter]);
			counter++;
			if (!all_blocks_polynomial(tempvars))
				return false;
		}
		for (const auto& vars : hom_variable_groups_)
			if (!all_blocks_polynomial(vars))
				return false;
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
		InvalidateDifferentiation();
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

		InvalidateDifferentiation();
		have_ordering_ = false;
		is_patched_ = false;
	}




	void System::AddHomVariableGroup(VariableGroup const& v)
	{
		hom_variable_groups_.push_back(v);
		InvalidateDifferentiation();
		have_ordering_ = false;
		is_patched_ = false;
		time_order_of_variable_groups_.push_back( VariableGroupType::Homogeneous);
	}





	void System::AddUngroupedVariable(Var const& v)
	{
		ungrouped_variables_.push_back(v);
		InvalidateDifferentiation();
		have_ordering_ = false;
		is_patched_ = false;
		time_order_of_variable_groups_.push_back( VariableGroupType::Ungrouped);
	}




	void System::AddUngroupedVariables(VariableGroup const& v)
	{
		ungrouped_variables_.insert( ungrouped_variables_.end(), v.begin(), v.end() );
		InvalidateDifferentiation();
		have_ordering_ = false;
		is_patched_ = false;
		time_order_of_variable_groups_.insert(time_order_of_variable_groups_.end(), v.size(), VariableGroupType::Ungrouped);
	}



 
	void System::AddImplicitParameter(Var const& v)
	{
		implicit_parameters_.push_back(v);
		InvalidateDifferentiation();
	}




	void System::AddImplicitParameters(VariableGroup const& v)
	{
		implicit_parameters_.insert( implicit_parameters_.end(), v.begin(), v.end() );
		InvalidateDifferentiation();
	}









	void System::AddParameter(Fn const& F)
	{
		explicit_parameters_.push_back(F);
		InvalidateDifferentiation();
	}



	void System::AddParameters(std::vector<Fn> const& v)
	{
		explicit_parameters_.insert( explicit_parameters_.end(), v.begin(), v.end() );
		InvalidateDifferentiation();
	}





	void System::AddSubfunction(Fn const& F)
	{
		PolyBlock().AddSubFunction(F);
		InvalidateDifferentiation();
	}



	void System::AddSubfunctions(std::vector<Fn> const& v)
	{
		for (auto const& f : v) PolyBlock().AddSubFunction(f);
		InvalidateDifferentiation();
	}





	void System::AddFunction(Fn const& F)
	{
		PolyBlock().AddFunction(F);
		InvalidateDifferentiation();
	}



	void System::AddFunction(Nd const& N, std::string const& name)
	{
		PolyBlock().AddFunction(Function::Make(N, name));
		InvalidateDifferentiation();
	}



	void System::AddFunctions(std::vector<Fn> const& v)
	{
		for (auto const& f : v) PolyBlock().AddFunction(f);
		InvalidateDifferentiation();
	}






	void System::AddConstant(Fn const& F)
	{
		PolyBlock().AddConstant(F);
		InvalidateDifferentiation();
	}


	void System::AddConstants(std::vector<Fn> const& v)
	{
		for (auto const& f : v) PolyBlock().AddConstant(f);
		InvalidateDifferentiation();
	}






	void System::AddPathVariable(Var const& v)
	{
		path_variable_ = v;
		InvalidateDifferentiation();
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
    	if (degs.empty())
    		return 0;   // a system with no functions has no degree bound
    	return *std::max_element(degs.begin(), degs.end());
    }


	std::vector<int> System::Degrees() const
	{
		std::vector<int> degs;
		for (auto const& b : blocks_)
		{
			auto d = std::visit([](auto const& blk){ return blk.Degrees(); }, b);
			degs.insert(degs.end(), d.begin(), d.end());
		}
		return degs;
	}


	std::vector<int> System::Degrees(VariableGroup const& vars) const
	{
		std::vector<int> degs;
		for (auto const& b : blocks_)
		{
			auto d = std::visit([&](auto const& blk){ return blk.Degrees(vars); }, b);
			degs.insert(degs.end(), d.begin(), d.end());
		}
		return degs;
	}


	//
	//  ExpandToFunctionTree -- build the pure function-tree twin of a block-composed system.
	//  A verification / interop oracle (see the header).  Scoped to the current block types.
	//
	namespace {

		// Build the function-tree node for a single linear form  sum_c M(r,c)*var_c + M(r,n),
		// where row r of M holds the (augmented) coefficients and `vars` are the ordered variable
		// nodes (column c <-> vars[c]); the last column is the constant / augmenting term.  Zero
		// coefficients are skipped to keep the tree compact (and exact: a skipped term is +0).
		Nd LinearFormNode(Mat<mpfr_complex> const& M, Eigen::Index r,
		                  VariableGroup const& vars, size_t num_vars)
		{
			Nd form = node::Float::Make(M(r, static_cast<Eigen::Index>(num_vars))); // constant term
			for (size_t c = 0; c < num_vars; ++c)
			{
				mpfr_complex const& coeff = M(r, static_cast<Eigen::Index>(c));
				if (coeff.real() == 0 && coeff.imag() == 0)
					continue;
				form = form + node::Float::Make(coeff) * vars[c];
			}
			return form;
		}

	} // anonymous namespace


	std::vector<Nd> System::NaturalFunctionsAsNodes() const
	{
		using namespace bertini::node;
		std::vector<Nd> out;

		auto const& vars = Variables(); // ordered variable nodes; block column c <-> vars[c]

		for (auto const& blk : blocks_)
		{
			std::visit([&](auto const& b)
			{
				using B = std::decay_t<decltype(b)>;

				if constexpr (std::is_same_v<B, blocks::PolynomialBlock>)
				{
					// already function-tree: take each Function's inner expression node.
					for (auto const& f : b.Functions())
						out.push_back(f->EntryNode());
				}
				else if constexpr (std::is_same_v<B, blocks::ProductsOfLinearsBlock>)
				{
					// f_i = prod_r ( row r of factor-matrix i . [vars ; 1] )
					const size_t n = b.NumVariables();
					if (static_cast<size_t>(vars.size()) != n)
						throw std::runtime_error("ExpandToFunctionTree: products-of-linears variable count mismatch");
					for (auto const& M : b.Factors())
					{
						Nd prod = nullptr;
						for (Eigen::Index r = 0; r < M.rows(); ++r)
						{
							Nd factor = LinearFormNode(M, r, vars, n);
							prod = prod ? (prod * factor) : factor;
						}
						out.push_back(prod ? prod : Nd(Integer::Make(1))); // empty product == 1
					}
				}
				else if constexpr (std::is_same_v<B, blocks::BlendBlock<System>>)
				{
					// H = sum_i c_i(t) * operand_i, each operand expanded to nodes (recursion).
					auto const& coeffs   = b.Coefficients();
					auto const& operands = b.Operands();
					const size_t k = b.NumFunctions();
					std::vector<Nd> blended(k, nullptr);
					for (size_t i = 0; i < operands.size(); ++i)
					{
						std::vector<Nd> fi = operands[i]->NaturalFunctionsAsNodes();
						if (fi.size() < k)
							throw std::runtime_error("ExpandToFunctionTree: blend operand has too few functions");
						for (size_t j = 0; j < k; ++j)
						{
							Nd term = coeffs[i] * fi[j];
							blended[j] = blended[j] ? (blended[j] + term) : term;
						}
					}
					for (auto& f : blended)
						out.push_back(f ? f : Nd(Integer::Make(0)));
				}
				else if constexpr (std::is_same_v<B, blocks::RandomizationBlock<System>>)
				{
					// g_i = sum_j c_ij * f_j * prod_g h_g^{(D_{i,g} - d_{j,g})}, the operand functions
					// expanded recursively and the homogenizing-variable powers folded back in (so the
					// expansion matches the block's homogenized evaluation).
					std::vector<Nd> fj = b.Operand()->NaturalFunctionsAsNodes();      // N nodes
					auto const& R   = b.RandomizationMatrix();                        // n x N
					auto const& tgt = b.TargetMultidegrees();
					auto const& omd = b.OperandMultidegrees();
					auto const& homvars = b.HomVars();
					const bool hom = b.IsHomogenized();
					const size_t n = b.NumFunctions();
					const size_t N = fj.size();
					for (size_t i = 0; i < n; ++i)
					{
						Nd gi = nullptr;
						for (size_t j = 0; j < N; ++j)
						{
							mpfr_complex const& c = R(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
							if (c.real() == 0 && c.imag() == 0)
								continue;
							Nd term = Float::Make(c) * fj[j];
							if (hom)
								for (size_t g = 0; g < homvars.size(); ++g)
								{
									const int e = tgt[i][g] - omd[j][g];
									if (e > 0)
										term = term * pow(homvars[g], e);
								}
							gi = gi ? (gi + term) : term;
						}
						out.push_back(gi ? gi : Nd(Integer::Make(0)));
					}
				}
				else if constexpr (std::is_same_v<B, blocks::LinearFormsBlock>)
				{
					// each row is one affine linear form  sum_c M(r,c)*vars[c] (+ constant).  Affine:
					// the trailing column is the constant.  Homogeneous (post-Homogenize): every column
					// is a variable column (the old constant is now the homogenizing-variable coeff).
					auto const& M = b.Coefficients();
					const size_t n = b.NumVariables();
					if (static_cast<size_t>(vars.size()) != n)
						throw std::runtime_error("ExpandToFunctionTree: linear-forms variable count mismatch");
					for (Eigen::Index r = 0; r < M.rows(); ++r)
					{
						if (b.IsHomogenized())
						{
							Nd form = nullptr;
							for (size_t c = 0; c < n; ++c)
							{
								mpfr_complex const& coeff = M(r, static_cast<Eigen::Index>(c));
								if (coeff.real() == 0 && coeff.imag() == 0)
									continue;
								Nd term = node::Float::Make(coeff) * vars[c];
								form = form ? (form + term) : term;
							}
							out.push_back(form ? form : Nd(Integer::Make(0)));
						}
						else
						{
							out.push_back(LinearFormNode(M, r, vars, n));   // augmented: last col is the constant
						}
					}
				}
				else // any future block
				{
					throw std::runtime_error("ExpandToFunctionTree: block type not yet supported");
				}
			}, blk);
		}

		return out;
	}


	System System::ExpandToFunctionTree() const
	{
		// expand THIS system's blocks to nodes first (reads the current block structure)...
		std::vector<Nd> nodes = NaturalFunctionsAsNodes();

		// ...then build the twin from a copy (same variables / groups / hom vars / path variable /
		// patch / ordering), with all blocks replaced by one PolynomialBlock of those nodes.
		System result = *this;
		result.ClearBlocks();
		for (auto const& f : nodes)
			result.AddFunction(f);
		result.InvalidateDifferentiation();
		return result;
	}


	//
	//  Randomize -- square up an overdetermined system (see the header).  Construction (degrees,
	//  sorting, the coefficient matrix) happens here, on a copy, where System is complete; the
	//  RandomizationBlock just stores the finished matrix and multidegrees and evaluates.
	//
	namespace {

		// operand->Degrees(group_g)[j] gathered into operand_multidegrees[j][g].
		std::vector<std::vector<int>> OperandMultidegrees(System const& operand)
		{
			auto groups = operand.VariableGroups();
			const size_t G = groups.size();
			const size_t N = operand.NumNaturalFunctions();
			std::vector<std::vector<int>> md(N, std::vector<int>(G, 0));
			for (size_t g = 0; g < G; ++g)
			{
				auto dg = operand.Degrees(groups[g]);            // length N: degree of each function in group g
				for (size_t j = 0; j < N && j < dg.size(); ++j)
					md[j][g] = dg[j];
			}
			return md;
		}

	} // anonymous namespace


	System System::AssembleRandomized(std::shared_ptr<System> operand, Mat<mpfr_complex> coefficients) const
	{
		const size_t G = operand->NumVariableGroups();
		const size_t N = operand->NumNaturalFunctions();
		const size_t n = static_cast<size_t>(coefficients.rows());

		if (static_cast<size_t>(coefficients.cols()) != N)
			throw std::runtime_error("Randomize: coefficient matrix column count must equal the number of natural functions.");

		auto operand_md = OperandMultidegrees(*operand);

		// Row i's target multidegree is, per group, the largest degree among the functions actually
		// combined into it (those with a nonzero coefficient).  This makes every h-power deficit
		// D_{i,g} - d_{j,g} >= 0, and -- with the descending sort the auto path uses -- equal to the
		// row's own leading-function degree, so the path count is minimal.
		std::vector<std::vector<int>> target_md(n, std::vector<int>(G, 0));
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < N; ++j)
			{
				mpfr_complex const& c = coefficients(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
				if (c.real() == 0 && c.imag() == 0)
					continue;
				for (size_t g = 0; g < G; ++g)
					target_md[i][g] = std::max(target_md[i][g], operand_md[j][g]);
			}

		blocks::RandomizationBlock<System> block(operand, std::move(coefficients),
		                                         std::move(target_md), std::move(operand_md), G);

		System result = *this;          // share variables / groups / path variable / ordering
		result.ClearBlocks();           // drop this system's own functions...
		result.AddBlock(std::move(block));  // ...the randomized rows come from the block
		result.InvalidateDifferentiation();
		return result;
	}


	System System::Randomize() const
	{
		const size_t G = NumVariableGroups();
		const size_t N = NumNaturalFunctions();
		const size_t n = NumVariables() - NumHomVariableGroups();

		if (N < n)
			throw std::runtime_error("Randomize: system is underdetermined (fewer functions than variables), so it has no isolated solutions to capture.");

		auto operand = std::make_shared<System>(*this);

		Mat<mpfr_complex> R(static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(N));

		if (G == 1)
		{
			// single affine group: sort the operand's functions by descending degree, then R = [I | C].
			// The identity block makes g_i carry f_i with coefficient 1 (degree d_i); the random tail
			// C folds the lower-degree functions in, padded by hom-var powers.  deg g_i = d_i, so the
			// total-degree path count is the product of the n largest degrees -- optimal.
			operand->ReorderFunctionsByDegreeDecreasing();
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < N; ++j)
				{
					if (j < n)
						R(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
							(i == j) ? mpfr_complex(1) : mpfr_complex(0);
					else
						R(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
							bertini::multiprecision::RandomComplex(DefaultPrecision());
				}
		}
		else
		{
			// several variable groups: multidegrees are only partially ordered, so use a dense random
			// R with a common (componentwise-max) target multidegree -- correct, and optimal when the
			// functions share a multidegree.
			for (Eigen::Index r = 0; r < R.rows(); ++r)
				for (Eigen::Index c = 0; c < R.cols(); ++c)
					R(r, c) = bertini::multiprecision::RandomComplex(DefaultPrecision());
		}

		return AssembleRandomized(operand, std::move(R));
	}


	System System::Randomize(Mat<mpfr_complex> const& R) const
	{
		if (static_cast<size_t>(R.cols()) != NumNaturalFunctions())
			throw std::runtime_error("Randomize: supplied matrix must have one column per natural function of the system.");
		auto operand = std::make_shared<System>(*this);  // functions kept in their current order
		return AssembleRandomized(operand, R);
	}


	Mat<mpfr_complex> System::RandomizationMatrix() const
	{
		for (auto const& b : blocks_)
			if (auto const* rb = std::get_if<blocks::RandomizationBlock<System>>(&b))
				return rb->RandomizationMatrix();
		throw std::runtime_error("RandomizationMatrix: this system has no randomization block (it was not produced by Randomize()).");
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
			re_ordered_functions[ind] = PolyBlock().Functions()[iter];
			ind++;
		}

		swap(PolyBlock().Functions(), re_ordered_functions);
		InvalidateDifferentiation();
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
			re_ordered_functions[ind] = PolyBlock().Functions()[iter];
			ind++;
		}

		swap(PolyBlock().Functions(), re_ordered_functions);
		InvalidateDifferentiation();
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

		InvalidateDifferentiation();
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
			InvalidateDifferentiation();
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
		if (auto* p = PolyBlockPtr())
			p->SimplifyFunctions();

		InvalidateDifferentiation();
	}



	void System::SimplifyDerivatives() const
	{
		SyncPolyBlock();
		if (auto* p = PolyBlockPtr())
		{
			if (!p->IsDifferentiated()) p->Differentiate();
			p->SimplifyDerivatives();
		}
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

	void System::Describe(std::ostream& out, bool verbose) const
	{
		// --- variables ---
		out << NumVariableGroups() << (NumVariableGroups() == 1 ? " variable group:\n" : " variable groups:\n");
		{
			auto counter = 0;
			for (const auto& grp : variable_groups_)
			{
				out << "  group " << counter++ << ": ";
				for (auto const& v : grp) out << *v << " ";
				out << "\n";
			}
		}
		if (!hom_variable_groups_.empty())
		{
			out << NumHomVariableGroups() << " projective variable groups:\n";
			auto counter = 0;
			for (const auto& grp : hom_variable_groups_)
			{
				out << "  group " << counter++ << ": ";
				for (auto const& v : grp) out << *v << " ";
				out << "\n";
			}
		}
		if (NumHomVariables() != 0)
		{
			out << "  homogenizing variables: ";
			for (const auto& v : homogenizing_variables_) out << *v << " ";
			out << "\n";
		}
		if (!ungrouped_variables_.empty())
		{
			out << "  ungrouped variables: ";
			for (const auto& v : ungrouped_variables_) out << *v << " ";
			out << "\n";
		}

		// --- functions, block by block ---
		out << "\n" << NumNaturalFunctions() << (NumNaturalFunctions() == 1 ? " function:\n" : " functions:\n");
		VariableGroup vars;
		try { vars = VariableOrdering(); } catch (...) { /* unordered/malformed: print without var names */ }
		size_t row = 0;
		for (auto const& blk : blocks_)
			std::visit([&](auto const& b){ b.Describe(out, row, vars, verbose); }, blk);

		// --- parameters / constants (only when present) ---
		if (NumParameters())
		{
			out << "\n" << NumParameters() << " explicit parameters:\n";
			for (const auto& p : explicit_parameters_)
				out << "  " << p->name() << " = " << p->EntryNode() << "\n";
		}
		if (NumConstants())
		{
			out << "\n" << NumConstants() << " constants:\n";
			for (const auto& c : PolyBlockPtr()->ConstantSubfunctions())
				out << "  " << c->name() << " = " << c->EntryNode() << "\n";
		}

		// --- path variable / patch (only the informative bits) ---
		if (path_variable_)
			out << "\npath variable: " << path_variable_->name() << "\n";
		if (IsPatched())
		{
			if (verbose)
				out << "\n" << patch_;
			else
				out << "\npatched (" << NumPatches() << (NumPatches() == 1 ? " patch)\n" : " patches)\n");
		}
	}


	std::ostream& operator<<(std::ostream& out, const bertini::System & s)
	{
		s.Describe(out, /*verbose=*/false);
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
		{
			auto& lhsf = PolyBlock().Functions();
			auto const& rhsf = rhs.PolyFunctions();
			for (size_t ii = 0; ii < lhsf.size(); ++ii)
				lhsf[ii] = node::Function::Make(rhsf[ii]->EntryNode() + lhsf[ii]->EntryNode(), lhsf[ii]->name());
		}

		InvalidateDifferentiation();
		return *this;
	}

	const System operator+(System lhs, System const& rhs)
	{
		return lhs+=rhs;
	}


	System& System::operator*=(std::shared_ptr<node::Node> const& N)
	{
		// new wrappers, not SetRoot — see comment in operator+= above.
		for (auto& f : PolyBlock().Functions())
		{
			f = node::Function::Make( N * f->EntryNode(), f->name());
		}
		InvalidateDifferentiation();
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


	System MakeHomotopy(System const& target, System const& start,
	                    std::string const& path_variable_name,
	                    std::shared_ptr<node::Node> const& gamma)
	{
		auto t = node::Variable::Make(path_variable_name);
		auto g = gamma ? gamma
		               : std::static_pointer_cast<node::Node>(node::Rational::Make(node::Rational::Rand()));

		System homotopy;
		if (start.HasStructuredBlocks() || target.HasStructuredBlocks())
		{
			// A block-backed system (e.g. a products-of-linears start, or a randomized target)
			// cannot be fused into a node-arithmetic homotopy: operator+ / operator* only combine
			// the PolynomialBlock functions and silently ignore structured blocks.  So whenever
			// EITHER side carries a structured block, combine the two systems with a blend block:
			// H = (1-t)*target + gamma*t*start, evaluated by blending whole Systems.  The homotopy
			// carries target's variable structure and patch; the blend contributes the natural
			// rows.  ClearBlocks drops the shell's own function blocks (a structured target's rows
			// live in a structured block, not a PolynomialBlock, so ClearFunctions would leave them
			// to be evaluated a second time alongside the blend).  Mirrors policy::CloneGiven::FormHomotopy.
			homotopy = target;
			homotopy.ClearBlocks();
			homotopy.AddPathVariable(t);
			std::vector<std::shared_ptr<node::Node>> coeffs{ 1 - t, g * t };
			std::vector<std::shared_ptr<const System>> operands{
				std::make_shared<System>(target),
				std::make_shared<System>(start) };
			homotopy.AddBlock(blocks::BlendBlock<System>(t, std::move(coeffs), std::move(operands)));
		}
		else
		{
			homotopy = (1-t)*target + g*t*start;
			homotopy.AddPathVariable(t);
		}
		return homotopy;
	}


	System MakeMovingHomotopy(System const& fixed, System const& start_moving, System const& end_moving,
	                          std::string const& path_variable_name,
	                          std::shared_ptr<node::Node> const& gamma)
	{
		if (start_moving.NumNaturalFunctions() != end_moving.NumNaturalFunctions())
			throw std::runtime_error("MakeMovingHomotopy: start_moving and end_moving must have the same number of functions (they are the two endpoints of the moving rows).");
		if (fixed.NumVariables() != start_moving.NumVariables() || fixed.NumVariables() != end_moving.NumVariables())
			throw std::runtime_error("MakeMovingHomotopy: fixed, start_moving and end_moving must share the same variable structure.");
		if (start_moving.HavePathVariable() || end_moving.HavePathVariable() || fixed.HavePathVariable())
			throw std::runtime_error("MakeMovingHomotopy: the fixed and moving systems must not already have a path variable.");

		auto t = node::Variable::Make(path_variable_name);
		auto g = gamma ? gamma
		               : std::static_pointer_cast<node::Node>(node::Rational::Make(node::Rational::Rand()));

		// Keep the fixed system's blocks as sibling blocks (do NOT clear them): they are autonomous,
		// so they are evaluated once per point and contribute nothing to dH/dt as the moving rows
		// slide.  Append a single blend block that moves only the moving rows:
		//   moving = (1-t)*end_moving + gamma*t*start_moving   (t=1 -> gamma*start, t=0 -> end).
		// Mirrors MakeHomotopy's blend branch, but blends only the moving operands instead of whole
		// systems, so the fixed equations are never duplicated or scaled.
		System homotopy = fixed;
		homotopy.AddPathVariable(t);
		std::vector<std::shared_ptr<node::Node>> coeffs{ 1 - t, g * t };
		std::vector<std::shared_ptr<const System>> operands{
			std::make_shared<System>(end_moving),
			std::make_shared<System>(start_moving) };
		homotopy.AddBlock(blocks::BlendBlock<System>(t, std::move(coeffs), std::move(operands)));
		return homotopy;
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
		if (auto* p = PolyBlockPtr()) p->Reset();
	}

	void System::ResetJacobian() const
	{
		if (auto* p = PolyBlockPtr()) p->Reset();
	}

	void System::ResetTimeDerivatives() const
	{
		if (auto* p = PolyBlockPtr()) p->Reset();
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
