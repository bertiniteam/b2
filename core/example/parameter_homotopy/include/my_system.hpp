
#pragma once

#include <bertini2/system.hpp>

namespace demo{


	using Node = std::shared_ptr<bertini::node::Node>;
	using Variable = std::shared_ptr<bertini::node::Variable>;

	using dbl = bertini::dbl;
	using mpfr = bertini::mpfr;

	auto MakeStep1Parameters()
	{
		using bertini::MakeVariable;

		// make symbolic objects for the parameters
		auto param_A = MakeFloat(bertini::RandomComplex(30));
		auto param_B = MakeFloat(bertini::RandomComplex(30));
		auto param_C = MakeFloat(bertini::RandomComplex(30));
		auto param_D = MakeFloat(bertini::RandomComplex(30));

		return std::vector<Node>{param_A, param_B, param_C, param_D};
	}


	template<typename ParamContT>
	auto MakeStep2Parameters(ParamContT const& step1_params, Node const& time)
	{

		using bertini::MakeVariable;

		std::vector<Node> steptwo_param_funcs;
		std::vector<Variable> steptwo_params;

		std::string suffix = "A";
		for (const auto& p: step1_params)
		{
			auto new_param = MakeVariable("param_"+suffix);
			steptwo_params.push_back(new_param);
			steptwo_param_funcs.push_back((1-time)*new_param + time*p);
			++suffix[0];
		}

		return std::make_tuple(steptwo_param_funcs, steptwo_params);
	}


	template <typename ParamContT>
	auto ConstructSystem(ParamContT const& params)
	{
		using bertini::MakeVariable;

	    auto x1 = MakeVariable("x1");
	    auto x2 = MakeVariable("x2");
	    auto x3 = MakeVariable("x3");
	    auto x4 = MakeVariable("x4");

	    auto f1 = x1*x1*x1*params[0] + x1*x1*x2*params[1] + x1*x2*x2*params[2] + x1*x3*x3*params[3] + x1*x4*x4*params[0]
	     + x1*params[1]+ x2*x2*x2*params[2] + x2*x3*x3*params[3] + x2*x4*x4*params[0] + x2*params[1] + 1;
	     
	    auto f2 = x1*x1*x1*params[2] + x1*x1*x2*params[3] + x1*x2*x2*params[0] + x1*x3*x3*params[1] + x1*x4*x4*params[2]
	     + x1*params[3] + x2*x2*x2*params[0] + x2*x3*x3*params[1] + x2*x4*x4*params[2] + x2*params[3] - 1;
	     
	    auto f3 = x1*x1*x3*params[0] + x1*x2*x3*params[1] + x2*x2*x3*params[2] + x3*x3*x3*params[3] + x3*x4*x4*params[0] + x3*params[1] + 2;

	    auto f4 = x1*x1*x4*params[2] + x1*x2*x4*params[3] + x2*x2*x4*params[0] + x3*x3*x4*params[1] + x4*x4*x4*params[2] + x4*params[3] - 3;

	    // make an empty system
	    bertini::System Sys;

	    // add the functions.  we could elide the `auto` construction above and construct directly into the system if we wanted
	    Sys.AddFunction(f1);
	    Sys.AddFunction(f2);
	    Sys.AddFunction(f3);
	    Sys.AddFunction(f4);

	    // make an affine variable group
	    bertini::VariableGroup vg{x1, x2, x3, x4};
	    Sys.AddVariableGroup(vg);

	    return Sys;
	}



	auto ConstructStart(bertini::System const& sys)
	{
	    return bertini::start_system::TotalDegree(sys);
	}

	template <typename StartT>
	auto ConstructHomotopy(bertini::System const& target_sys, StartT const& start_sys)
	{
		using bertini::MakeVariable;

		auto t = MakeVariable("t");

		auto gamma = bertini::MakeRational(bertini::node::Rational::Rand());

		return (1-t)*target_sys + gamma*t*start_sys;
	}


} // namespace demo

