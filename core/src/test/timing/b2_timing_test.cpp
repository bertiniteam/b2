#include "bertini.hpp"

#include <boost/timer/timer.hpp>


using System = bertini::System;
using Var = std::shared_ptr<bertini::Variable>;
using VariableGroup = bertini::VariableGroup;


template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;



void simple_single_variable();
void complicated_single_variable();

int main(int argc, char** argv)
{	
	simple_single_variable();

	complicated_single_variable();

	return 0;
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

	boost::timer::auto_cpu_timer timing_guy;

	auto J = S.Jacobian(v,dbl(0));
	auto f_values = S.Eval(v,dbl(0));
	for (unsigned ii=0; ii < 1000000; ii++)
	{
		f_values = S.Eval(v,dbl(0));
		// J = S.Jacobian<dbl>();
		// std::cout << J << "\n";
	}

	std::cout << S << "\n";
}


void complicated_single_variable()
{
	boost::timer::auto_cpu_timer t;

	std::cout << "complicated_single_variable" << "\n";
}


