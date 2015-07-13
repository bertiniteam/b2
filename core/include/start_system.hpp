//This file is part of Bertini 2.0.
//
//start_system.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//start_system.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with start_system.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  start_system.hpp
//
//  copyright 2015
//  Daniel Brake
//  University of Notre Dame
//  ACMS
//  Spring, Summer 2015


#include "system.hpp"



namespace bertini 
{
	namespace start_system{

		
		template<typename T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
		template<typename T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
			

		class StartSystem : public System
		{

		public:
			

			virtual size_t NumStartPoints() const = 0;

			template<typename T>
			Vec<T> StartPoint(size_t index) const
			{
				return GenerateStartPoint(T(),index);
			}

		private:
			virtual Vec<dbl> GenerateStartPoint(dbl,size_t index) const = 0;
			virtual Vec<mpfr> GenerateStartPoint(mpfr,size_t index) const = 0;

		};


		class TotalDegree : public StartSystem
		{
		public:
			TotalDegree() = default;
			virtual ~TotalDegree() = default;
			
			/**
			 Constructor for making a total degree start system from a polynomial system
			*/
			TotalDegree(System const& s);


			mpfr RandomValue(size_t index)
			{
				return random_values_[index]->Eval<mpfr>();
			}


			std::vector<std::shared_ptr<Rational> > const& RandomValues()
			{
				return random_values_;
			}

			size_t NumStartPoints() const override;


		private:


			Vec<dbl> GenerateStartPoint(dbl,size_t index) const override;
			Vec<mpfr> GenerateStartPoint(mpfr,size_t index) const override;

			std::vector<std::shared_ptr<Rational> > random_values_; ///< stores the random values for the start functions.  x^d-r, where r is stored in this vector.
			std::vector<int> degrees_;
		};
	}
}
