//This file is part of Bertini 2.
//
//mhom.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mhom.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with mhom.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// daniel brake, university of notre dame

#include "bertini2/system/start/mhom.hpp"


BOOST_CLASS_EXPORT(bertini::start_system::MHomogeneous);


namespace bertini {

	namespace start_system {

		// constructor for MHomogeneous start system, from any other *suitable* system.
		MHomogeneous::MHomogeneous(System const& s)
		{

			// if (s.NumHomVariableGroups() > 0)
			// 	throw std::runtime_error("a homogeneous variable group is present.  currently unallowed");


			// if (s.NumTotalFunctions() != s.NumVariables())
			// 	throw std::runtime_error("attempting to construct total degree start system from non-square target system");

			// if (s.HavePathVariable())
			// 	throw std::runtime_error("attempting to construct total degree start system, but target system has path varible declared already");			

			// if (!s.IsPolynomial())
			// 	throw std::runtime_error("attempting to construct total degree start system from non-polynomial target system");

			CreateDegreeMatrix(s);

			GenerateValidPartitions(s);


			// auto deg = s.Degrees();
			// for (auto iter : deg)
			// 	degrees_.push_back((size_t) iter);

			// CopyVariableStructure(s);

			// random_values_.resize(s.NumFunctions());


			// for (unsigned ii = 0; ii < s.NumFunctions(); ++ii)
			// 	random_values_[ii] = MakeRational(node::Rational::Rand());

			// // by hypothesis, the system has a single variable group.
			// VariableGroup v = this->AffineVariableGroup(0);

			// for (auto iter = v.begin(); iter!=v.end(); iter++)
			// 	AddFunction(pow(*iter,(int) *(degrees_.begin() + (iter-v.begin()))) - random_values_[iter-v.begin()]); 
			
			// if (s.IsHomogeneous())
			// 	Homogenize();

			// if (s.IsPatched())
			// 	CopyPatches(s);

		}// total degree constructor

		
		// MHomogeneous& MHomogeneous::operator*=(Nd const& n)
		// {
		// 	*this *= n;
		// 	return *this;
		// }
		
		

		// unsigned long long MHomogeneous::NumStartPoints() const
		// {
		// 	unsigned long long num_start_points = 1;
		// 	for (auto iter : degrees_)
		// 		num_start_points*=iter;
		// 	return num_start_points;
		// }

		void MHomogeneous::CreateDegreeMatrix(System const& s)
		{
			degree_matrix_ = Mat<int>(s.NumFunctions(),s.NumHomVariableGroups());
			auto var_groups = s.HomVariableGroups();
			int col_count = 0;
			for(std::vector<VariableGroup>::iterator it = var_groups.begin(); it != var_groups.end(); ++it) 
			{
  				std::vector<int> degs = s.Degrees(*it);

  				for(int ii = 0; ii <= degs.size() - 1; ++ii)
  				{
  					degree_matrix_(ii,col_count) = degs[ii];
  				}
  				col_count++;
			}

		}

		void MHomogeneous::GenerateValidPartitions(System const& s)
		{
		// 	int row = 0;
		// 	//  s.NumVariables();//  m=PPD->num_hom_var_gp+PPD->num_var_gp;  //m=# var gps; we cycle through row and col in degree table mhomDeg.
		// 	int count=0;
		// 	//int old_current_part_row = -1;
		// 	//int bad_choice = 0;
		// 	std::vector<int> current_partition;  
		// 	//int *current_part = (int *)bmalloc(PPD->num_funcs * sizeof(int));  //holds our current choice of partition
		// 	int *var_gp_ctr = (int *)bmalloc(s.NumVariables() * sizeof(int));
		// 	// int *var_gp_ctr = (int *)bmalloc(m * sizeof(int));  //tracks # of each var gp type chosen so far

		// 	auto size_of_each_var_gp = s.VariableGroupSizes();
		// 	for(int ii = 0; ii < s.NumVariables(); ++ii)
		// 	{
		// 		var_gp_ctr[ii] = size_of_each_var_gp[ii]
		// 	}			    

		// 	  for (i = 0; i < PPD->num_funcs; i++)  //for each function
		// 	  {
		// 	      current_partition.push_back(-1);
		// 	  }

			    
		// 	  while (row > -1)  // Algorithm will move up and down rows, kicking out to row=-1 at end
		// 	  {
		// 	    old_current_part_row = current_part[row];  //Hang on to previous choice of column for this row, in case we are done with this row.
		// 	    current_part[row] = choose_col_in_row(mhomDeg,var_gp_ctr,row,current_part[row],m);  //Pick next column (var gp) for the current row (func)
			      
		// 	    if (current_part[row] == m) // means we have exhausted all good columns for the current row, so we go back up a row
		// 	    {
		// 	      //      var_gp_ctr[old_current_part_row] = var_gp_ctr[old_current_part_row]+1;  //done with this row, so increment counter from previously chosen column for this row
		// 	      row = row-1;  //go back up a row
		// 	      bad_choice = 1;
		// 	    }
		// 	    else  //found a good choice of column for this row!
		// 	    {
		// 	      row = row+1;  //move on to next row!
		// 	      if (row < PPD->num_funcs)
		// 	        current_part[row] = -1; //since we are starting a new row, we start with the left-most entry (choose_col_in_row() first increments col)
		// 	    }
			     
		// 	    if ((row == PPD->num_funcs) && (!bad_choice))  //We have reached the final row with a good partition!
		// 	    {
		// 	      // NEEDED   count += generateFromPartition_d(PPD, coeff, patchCoeff, current_part, P, mhomDeg, OUT);  //This line generates all startpoints for this partition
		// 	     // var_gp_ctr[current_part[row-1]] = var_gp_ctr[current_part[row-1]] + 1;  //done with this row, so increment counter from previously chosen column for this row
		// 	      row = row - 1; //put the counter back on the last row to try to move to the next column
		// 	    }
		// 	    bad_choice=0;
		// 	  }
		// 	  return count;
		// 	}

		}


		
		// Vec<dbl> MHomogeneous::GenerateStartPoint(dbl,unsigned long long index) const
		// {
		// 	Vec<dbl> start_point(NumVariables());
		// 	auto indices = IndexToSubscript(index, degrees_);

		// 	unsigned offset = 0;
		// 	if (IsPatched())
		// 	{
		// 		start_point(0) = dbl(1);
		// 		offset = 1;
		// 	}

		// 	for (size_t ii = 0; ii< NumNaturalVariables(); ++ii)
		// 		start_point(ii+offset) = exp( std::acos(-1) * dbl(0,2) * double(indices[ii]) / double(degrees_[ii])  ) * pow(random_values_[ii]->Eval<dbl>(), double(1) / double(degrees_[ii]));

		// 	if (IsPatched())
		// 		RescalePointToFitPatchInPlace(start_point);

		// 	return start_point;
		// }


		// Vec<mpfr> MHomogeneous::GenerateStartPoint(mpfr,unsigned long long index) const
		// {
		// 	Vec<mpfr> start_point(NumVariables());
		// 	auto indices = IndexToSubscript(index, degrees_);

		// 	unsigned offset = 0;
		// 	if (IsPatched())
		// 	{
		// 		start_point(0) = mpfr(1);
		// 		offset = 1;
		// 	}

		// 	auto two_i = mpfr(0,2);
		// 	auto one = mpfr(1);

		// 	for (size_t ii = 0; ii< NumNaturalVariables(); ++ii)
		// 		start_point(ii+offset) = exp( acos( mpfr_float(-1) ) * two_i * mpfr_float(indices[ii]) / mpfr_float(degrees_[ii])  ) * pow(random_values_[ii]->Eval<mpfr>(), one / degrees_[ii]);

		// 	if (IsPatched())
		// 		RescalePointToFitPatchInPlace(start_point);


		// 	return start_point;
		// }

		// inline
		// MHomogeneous operator*(MHomogeneous td, std::shared_ptr<node::Node> const& n)
		// {
		// 	td *= n;
		// 	return td;
		// }

	} // namespace start_system
} //namespace bertini
