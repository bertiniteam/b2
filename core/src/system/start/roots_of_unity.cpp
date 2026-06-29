//This file is part of Bertini 2.
//
//roots_of_unity.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//roots_of_unity.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with roots_of_unity.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

#include "bertini2/system/start/roots_of_unity.hpp"

#include <boost/math/constants/constants.hpp>


BOOST_CLASS_EXPORT(bertini::start_system::RootsOfUnity);


namespace bertini {
	using namespace bertini::node;

	namespace start_system {

		// constructor for RootsOfUnity start system, from any other *suitable* system.
		RootsOfUnity::RootsOfUnity(System const& s)
		{
			SanityChecks(s);
			CopyDegrees(s);
			CopyVariableStructure(s);
			SeedRandomValues(static_cast<int>(s.NumNaturalFunctions()));
			GenerateFunctions();

			if (s.IsHomogeneous())
				Homogenize();

			if (s.IsPatched())
				CopyPatches(s);
		}// roots of unity constructor


		RootsOfUnity& RootsOfUnity::operator*=(Nd const& n)
		{
			System::operator*=(n);
			return *this;
		}



		unsigned long long RootsOfUnity::NumStartPoints() const
		{
			unsigned long long num_start_points = 1;
			for (const auto& iter : degrees_)
				num_start_points*=iter;
			return num_start_points;
		}



		Vec<complex_dbl> RootsOfUnity::GenerateStartPoint(complex_dbl,unsigned long long index) const
		{
			Vec<complex_dbl> start_point(NumVariables());
			auto indices = IndexToSubscript(index, degrees_);

			unsigned offset = 0;
			if (IsPatched())
			{
				start_point(0) = complex_dbl(1);
				offset = 1;
			}

			// authoritative pi (see issue #156 / special_number.cpp), not acos(-1)
			auto two_i_pi = boost::math::constants::pi<double>() * complex_dbl(0,2);

			for (size_t ii = 0; ii< NumNaturalVariables(); ++ii)
				start_point(static_cast<Eigen::Index>(ii+offset)) = exp( two_i_pi * static_cast<double>(indices[ii]) / static_cast<double>(degrees_[ii])  ) * pow(random_values_[ii]->Value<complex_dbl>(), 1.0 / static_cast<double>(degrees_[ii]));

			if (IsPatched())
				RescalePointToFitPatchInPlace(start_point);

			return start_point;
		}


		Vec<complex_mp> RootsOfUnity::GenerateStartPoint(complex_mp,unsigned long long index) const
		{
			using bertini::ThreadPrecision;

			Vec<complex_mp> start_point(NumVariables()); // make the value we're returning
			auto indices = IndexToSubscript(index, degrees_); // get the position of it -- used in the angle of the coordinates of the produced point.

			unsigned offset = 0;
			if (IsPatched())
			{
				start_point(0) = complex_mp(1,0,ThreadPrecision());
				offset = 1;
			}

// TODO: this code should be cleaned up after issue 308 is solved -- namely, the two precision adjustment calls should be removed.  They're only necessary because prec16 / ulonglog = prec19.

			auto one = real_mp(1);
			// authoritative pi (see issue #156 / special_number.cpp), not acos(-1)
			complex_mp two_i_pi = complex_mp(0,2) * boost::math::constants::pi<real_mp>();
			for (size_t ii = 0; ii< NumNaturalVariables(); ++ii)
			{
				complex_mp a = exp( (two_i_pi * indices[ii]) / degrees_[ii]);
				complex_mp b = pow(random_values_[ii]->Value<complex_mp>(), one / degrees_[ii]);

				Precision(a,ThreadPrecision());
				Precision(b,ThreadPrecision());

				start_point(static_cast<Eigen::Index>(ii+offset)) = a*b;
			}

			if (IsPatched())
				RescalePointToFitPatchInPlace(start_point);

			return start_point;
		}

		inline
		RootsOfUnity operator*(RootsOfUnity td, std::shared_ptr<node::Node> const& n)
		{
			td *= n;
			return td;
		}

		void RootsOfUnity::SanityChecks(System const& s)
		{
			if (s.NumHomVariableGroups() > 0)
				throw std::runtime_error("a homogeneous variable group is present.  currently unallowed");

			if (s.NumTotalFunctions() != s.NumVariables())
				throw std::runtime_error("attempting to construct roots-of-unity start system from non-square target system");

			if (s.HavePathVariable())
				throw std::runtime_error("attempting to construct roots-of-unity start system, but target system has path varible declared already");

			if (s.NumVariableGroups() != 1)
				throw std::runtime_error("more than one affine variable group.  currently unallowed");

			if (!s.IsPolynomial())
				throw std::runtime_error("attempting to construct roots-of-unity start system from non-polynomial target system");
		}

		void RootsOfUnity::CopyDegrees(System const& s)
		{
			auto deg = s.Degrees();
			for (const auto& d : deg)
				degrees_.push_back(static_cast<size_t>(d));
		}


		void RootsOfUnity::SeedRandomValues(int num_functions)
		{
			random_values_.resize(static_cast<size_t>(num_functions));
			for (int ii = 0; ii < num_functions; ++ii)
				random_values_[static_cast<size_t>(ii)] = Rational::Make(node::Rational::Rand());
		}

		void RootsOfUnity::GenerateFunctions()
		{
			// by hypothesis, the system has a single variable group.
			auto v = this->AffineVariableGroup(0);
			for (auto iter = v.begin(); iter!=v.end(); iter++)
				AddFunction(pow(*iter,(int) *(degrees_.begin() + (iter-v.begin()))) - random_values_[static_cast<size_t>(iter-v.begin())]);
		}
	} // namespace start_system
} //namespace bertini
