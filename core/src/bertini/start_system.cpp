//This file is part of Bertini 2.0.
//
//start_system.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//start_system.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with start_system.cpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  start_system.cpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015


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


			auto deg = s.Degrees();
			for (auto iter : deg)
				degrees_.push_back((size_t)iter);

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
				AddFunction(pow(*iter,(int) *(degrees_.begin() + (iter-original_variables.begin()))) - random_values_[iter-original_variables.begin()]); ///< TODO: make this 1 be a random complex number.
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
			auto indices = IndexToSubscript(index, degrees_);
			for (size_t ii = 0; ii< NumVariables(); ++ii)
				start_point(ii) = exp( acos(-1.0) * dbl(0,2.0) * double(indices[ii]) / double(degrees_[ii])  ) * pow(random_values_[ii]->Eval<dbl>(), 1.0 / degrees_[ii]);
			return start_point;
		}


		Vec<mpfr> TotalDegree::GenerateStartPoint(mpfr,size_t index) const
		{
			Vec<mpfr> start_point(NumVariables());
			auto indices = IndexToSubscript(index, degrees_);
			for (size_t ii = 0; ii< NumVariables(); ++ii)
				start_point(ii) = exp( acos( mpfr_float(-1.0) ) * mpfr(0,2.0) * mpfr_float(indices[ii]) / mpfr_float(degrees_[ii])  ) * pow(random_values_[ii]->Eval<mpfr>(), mpfr_float(1.0) / degrees_[ii]);
			return start_point;
		}


	} // namespace start_system
} //namespace bertini
