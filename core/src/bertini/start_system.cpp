#include "start_system.hpp"



namespace bertini {


	System TotalDegreeStartSystem(System const& s)
	{

		

		if (s.NumFunctions() != s.NumVariables())
			throw std::runtime_error("attempting to construct total degree start system from non-square target system");

		if (s.HavePathVariable())
			throw std::runtime_error("target system has path varible declared already");

		if (s.NumVariableGroups() != 1)
			throw std::runtime_error("more than one variable group.  currently unallowed");

		if (s.NumHomVariableGroups() > 0)
			throw std::runtime_error("a homogeneous variable group is present.  currently unallowed");

		if (!s.IsPolynomial())
			throw std::runtime_error("attempting to construct total degree start system from non-polynomial target system");


		auto degrees = s.Degrees();

		auto original_varible_groups = s.variableGroups();

		auto original_variables = s.variables();



		System TD;

		TD.CopyVariableStructure(s);


		for (auto iter = original_variables.begin(); iter!=original_variables.end(); iter++)
		{
			TD.AddFunction(pow(*iter,*(degrees.begin() + (iter-original_variables.begin()))) - 1); ///< TODO: make this 1 be a random complex number.
		} 

		return TD;
	}
	
}
