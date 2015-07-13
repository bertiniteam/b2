#include "start_system.hpp"



namespace bertini {

	namespace start_system {

		// constructor for TotalDegree start system, from any other *suitable* system.
		TotalDegree::TotalDegree(System const& s)
		{
			if (s.NumFunctions() != s.NumVariables())
				throw std::runtime_error("attempting to construct total degree start system from non-square target system");

			if (s.HavePathVariable())
				throw std::runtime_error("attempting to construct total degree start system, but target system has path varible declared already");

			if (s.NumVariableGroups() != 1)
				throw std::runtime_error("more than one variable group.  currently unallowed");

			if (s.NumHomVariableGroups() > 0)
				throw std::runtime_error("a homogeneous variable group is present.  currently unallowed");

			if (!s.IsPolynomial())
				throw std::runtime_error("attempting to construct total degree start system from non-polynomial target system");


			degrees_ = s.Degrees();

			auto original_varible_groups = s.variableGroups();

			auto original_variables = s.Variables();


			random_values_.resize(s.NumFunctions());

			// auto prev_prec = boost::multiprecision::mpfr_float::default_precision();
			// boost::multiprecision::mpfr_float::default_precision(4000);
			for (unsigned ii = 0; ii < s.NumFunctions(); ++ii)
			{
				random_values_[ii] = std::make_shared<Rational>(bertini::Rational::Rand());
			}


			CopyVariableStructure(s);


			for (auto iter = original_variables.begin(); iter!=original_variables.end(); iter++)
			{
				AddFunction(pow(*iter,*(degrees_.begin() + (iter-original_variables.begin()))) - random_values_[iter-original_variables.begin()]); ///< TODO: make this 1 be a random complex number.
			} 

			// boost::multiprecision::mpfr_float::default_precision(prev_prec);

		}// total degree constructor

		
		
		

		size_t TotalDegree::NumStartPoints() const
		{
			size_t num_start_points = 1;
			for (auto iter : degrees_)
				num_start_points*=iter;
			return num_start_points;
		}


		
		Vec<dbl> TotalDegree::GenerateStartPoint(dbl,size_t index) const
		{
			Vec<dbl> start_point(NumVariables());

			return start_point;
		}


		Vec<mpfr> TotalDegree::GenerateStartPoint(mpfr,size_t index) const
		{
			Vec<mpfr> start_point(NumVariables());

			return start_point;
		}


	} // namespace start_system
} //namespace bertini
