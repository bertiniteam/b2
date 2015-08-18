#include "bertini.hpp"

#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <fstream>

#include <sstream>



using System = bertini::System;
using Var = std::shared_ptr<bertini::Variable>;
using VariableGroup = bertini::VariableGroup;


template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
void arbitrary(boost::filesystem::path const& file, unsigned num_iterations = 1000000);

void simple_single_variable();

int main(int argc, char** argv)
{	
	switch (argc)
	{
		case 1:
		{
			simple_single_variable();
			break;
		}
		case 2:
		{	
			boost::filesystem::path file(argv[1]);
			arbitrary<dbl>(file);
			break;
		}
		default:
		{
			std::cout << "timing testing not implemented with the number of arguments you gave\n";
			break;
		}
	}

	

	return 0;
}


template<typename T>
void arbitrary(boost::filesystem::path const& file, unsigned num_iterations)
{
	std::ifstream fin(file.c_str());

	std::stringstream buffer;
	buffer << fin.rdbuf();

	System S(buffer.str());

	Vec<T> variable_values(S.NumVariables());

	for (unsigned i = 0; i < variable_values.size(); ++i)
		variable_values(i) = bertini::rand_complex();


	

	if (S.HavePathVariable())
	{
		T time;
		time = bertini::rand_complex();

		auto J = S.Jacobian(variable_values,time);
		auto f_values = S.Eval(variable_values,time);

		std::cout << S << "\n";

		boost::timer::auto_cpu_timer timing_guy;
		for (unsigned ii=0; ii<num_iterations; ii++)
		{
			f_values = S.Eval(variable_values,time);
		}
	}
	else
	{
		auto J = S.Jacobian(variable_values);
		auto f_values = S.Eval(variable_values);

		std::cout << S << "\n";

		boost::timer::auto_cpu_timer timing_guy;
		for (unsigned ii=0; ii<num_iterations; ii++)
		{
			f_values = S.Eval(variable_values);
		}
	}


}

void simple_single_variable()
{
	

	Var x = std::make_shared<bertini::Variable>("x");
	Var t = std::make_shared<bertini::Variable>("t");


	Vec<dbl> v(1);
	v << dbl(1.614234164678871592346182,-0.48717825396548971237946);

	System S;
	S.AddUngroupedVariable(x);
	S.AddPathVariable(t);
	S.AddFunction(pow(x,2) + x*t + pow(t,2));

	std::cout << S << "\n";
	
	boost::timer::auto_cpu_timer timing_guy;

	auto J = S.Jacobian(v,dbl(0));
	auto f_values = S.Eval(v,dbl(0));
	for (unsigned ii=0; ii < 1000000; ii++)
	{
		f_values = S.Eval(v,dbl(0));
		// J = S.Jacobian<dbl>();
		// std::cout << J << "\n";
	}

	
}



