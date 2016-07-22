//This file is part of Bertini 2.
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
// Copyright(C) 2015, 2016 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame

#include "start_system.hpp"


BOOST_CLASS_EXPORT(bertini::start_system::TotalDegree);


namespace bertini {

	namespace start_system {

		// constructor for TotalDegree start system, from any other *suitable* system.
		TotalDegree::TotalDegree(System const& s)
		{

			if (s.NumHomVariableGroups() > 0)
				throw std::runtime_error("a homogeneous variable group is present.  currently unallowed");


			if (s.NumTotalFunctions() != s.NumVariables())
				throw std::runtime_error("attempting to construct total degree start system from non-square target system");

			if (s.HavePathVariable())
				throw std::runtime_error("attempting to construct total degree start system, but target system has path varible declared already");

			if (s.NumVariableGroups() != 1)
				throw std::runtime_error("more than one affine variable group.  currently unallowed");

			

			if (!s.IsPolynomial())
				throw std::runtime_error("attempting to construct total degree start system from non-polynomial target system");

			

			auto deg = s.Degrees();
			for (auto iter : deg)
				degrees_.push_back((size_t) iter);

			CopyVariableStructure(s);

			random_values_.resize(s.NumFunctions());

			// auto prev_prec = boost::multiprecision::DefaultPrecision();
			// boost::multiprecision::DefaultPrecision(4000);
			for (unsigned ii = 0; ii < s.NumFunctions(); ++ii)
				random_values_[ii] = std::make_shared<node::Rational>(node::Rational::Rand());

			// by hypothesis, the system has a single variable group.
			VariableGroup v = this->AffineVariableGroup(0);

			for (auto iter = v.begin(); iter!=v.end(); iter++)
				AddFunction(pow(*iter,(int) *(degrees_.begin() + (iter-v.begin()))) - random_values_[iter-v.begin()]); 
			if (s.IsHomogeneous())
				Homogenize();

			if (s.IsPatched())
				CopyPatches(s);

			gamma_ = std::make_shared<node::Rational>(node::Rational::Rand());
			// *this *= static_cast<std::shared_ptr<node::Node> >(gamma_);
			// boost::multiprecision::DefaultPrecision(prev_prec);
		}// total degree constructor

		
		TotalDegree& TotalDegree::operator*=(Nd const& n)
		{
			*this *= n;
			return *this;
		}
		
		

		mpz_int TotalDegree::NumStartPoints() const
		{
			mpz_int num_start_points = 1;
			for (auto iter : degrees_)
				num_start_points*=iter;
			return num_start_points;
		}


		
		Vec<dbl> TotalDegree::GenerateStartPoint(dbl,mpz_int index) const
		{
			Vec<dbl> start_point(NumVariables());
			auto indices = IndexToSubscript(index, degrees_);

			unsigned offset = 0;
			if (IsPatched())
			{
				start_point(0) = dbl(1);
				offset = 1;
			}

			for (size_t ii = 0; ii< NumNaturalVariables(); ++ii)
				start_point(ii+offset) = exp( std::acos(-1) * dbl(0,2) * double(indices[ii]) / double(degrees_[ii])  ) * pow(random_values_[ii]->Eval<dbl>(), double(1) / double(degrees_[ii]));

			if (IsPatched())
				RescalePointToFitPatchInPlace(start_point);

			return start_point;
		}


		Vec<mpfr> TotalDegree::GenerateStartPoint(mpfr,mpz_int index) const
		{
			Vec<mpfr> start_point(NumVariables());
			auto indices = IndexToSubscript(index, degrees_);

			unsigned offset = 0;
			if (IsPatched())
			{
				start_point(0) = mpfr(1);
				offset = 1;
			}

			auto two_i = mpfr(0,2);
			auto one = mpfr(1);

			for (size_t ii = 0; ii< NumNaturalVariables(); ++ii)
				start_point(ii+offset) = exp( acos( mpfr_float(-1) ) * two_i * mpfr_float(indices[ii]) / mpfr_float(degrees_[ii])  ) * pow(random_values_[ii]->Eval<mpfr>(), one / degrees_[ii]);

			if (IsPatched())
				RescalePointToFitPatchInPlace(start_point);


			return start_point;
		}

		inline
		TotalDegree operator*(TotalDegree td, std::shared_ptr<node::Node> const& n)
		{
			td *= n;
			return td;
		}

	} // namespace start_system
} //namespace bertini
