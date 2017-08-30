
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
		auto param_A = MakeVariable("param_A");
		auto param_B = MakeVariable("param_B");
		auto param_C = MakeVariable("param_C");
		auto param_D = MakeVariable("param_D");

		auto val_A = bertini::RandomMp();
		auto val_B = bertini::RandomMp();
		auto val_C = bertini::RandomMp();
		auto val_D = bertini::RandomMp();


		param_A->set_current_value<dbl>(val_A.convert_to<double>());
		param_B->set_current_value<dbl>(val_B.convert_to<double>());
		param_C->set_current_value<dbl>(val_C.convert_to<double>());
		param_D->set_current_value<dbl>(val_D.convert_to<double>());

		param_A->set_current_value<mpfr>(val_A);
		param_B->set_current_value<mpfr>(val_B);
		param_C->set_current_value<mpfr>(val_C);
		param_D->set_current_value<mpfr>(val_D);

		return std::vector<Node>{param_A, param_B, param_C, param_D};
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
	     + x1*params[1]+ x2*x2*x2*params[2] + x2*x3*x3*params[3] + x2*x4*x4*params[0] + x2*params[1];
	     
	    auto f2 = x1*x1*x1*params[2] + x1*x1*x2*params[3] + x1*x2*x2*params[0] + x1*x3*x3*params[1] + x1*x4*x4*params[2]
	     + x1*params[3] + x2*x2*x2*params[0] + x2*x3*x3*params[1] + x2*x4*x4*params[2] + x2*params[3];
	     
	    auto f3 = x1*x1*x3*params[0] + x1*x2*x3*params[1] + x2*x2*x3*params[2] + x3*x3*x3*params[3] + x3*x4*x4*params[0] + x3*params[1];

	    auto f4 = x1*x1*x4*params[2] + x1*x2*x4*params[3] + x2*x2*x4*params[0] + x3*x3*x4*params[1] + x4*x4*x4*params[2] + x4*params[3];

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

