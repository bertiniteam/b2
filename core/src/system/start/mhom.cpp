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
// Copyright(C) 2015 - 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// karleigh cameron, colorado state university
// Tim Hodges, Colorado State University
// silviana amethyst, university of wisconsin eau claire

#include "bertini2/system/start/mhom.hpp"


BOOST_CLASS_EXPORT(bertini::start_system::MHomogeneous);


namespace bertini 
{

	namespace start_system 
	{

		// constructor for MHomogeneous start system, from any other *suitable* system. System can be homogeneous or non-homogeneous. 
		MHomogeneous::MHomogeneous(System const& s)
		{

			if (s.NumTotalFunctions() != s.NumVariables())
				throw std::runtime_error("attempting to construct multi homogeneous start system from non-square target system");

			if (s.HavePathVariable())
				throw std::runtime_error("attempting to construct multi homogeneous start system, but target system has path varible declared already");			

			if (!s.IsPolynomial())
				throw std::runtime_error("attempting to construct multi homogeneous start system from non-polynomial target system");


			CreateDegreeMatrix(s);


			GenerateValidPartitions(s);
			
			CopyVariableStructure(s);
			
			
			linprod_matrix_ = Mat<std::shared_ptr<node::LinearProduct>>(degree_matrix_.rows(), degree_matrix_.cols());
			std::shared_ptr<node::Node> func;
			for (int ii = 0; ii < degree_matrix_.rows(); ++ii)
			{
				func = bertini::MakeInteger(1);
				
				for (int jj = 0; jj < s.NumHomVariableGroups(); ++jj)
				{
					if(degree_matrix_(ii,jj) != 0)
					{
						// Fill the linear product matrix
						linprod_matrix_(ii,jj) = std::make_shared<node::LinearProduct>(var_groups_[jj], degree_matrix_(ii,jj), true);
						
						func *= linprod_matrix_(ii,jj);
					}
				}
				for (int jj = s.NumHomVariableGroups(); jj < degree_matrix_.cols(); ++jj)
				{
					if(degree_matrix_(ii,jj) != 0)
					{
						// Fill the linear product matrix
						linprod_matrix_(ii,jj) = std::make_shared<node::LinearProduct>(var_groups_[jj], degree_matrix_(ii,jj));
						
						func *= linprod_matrix_(ii,jj);
					}
				}
				
				AddFunction(func);
			}


			 if (s.IsHomogeneous())
			 	Homogenize();

			 if (s.IsPatched())
			 	CopyPatches(s);

		}// M-Hom constructor

		
		MHomogeneous& MHomogeneous::operator*=(Nd const& n)
		{
			*this *= n;
			return *this;
		}
		
		
		/**
			\brief Function to calculate the total number of start points. 

			## Input: 

			## Output:
					unsigned long long : Number representing the total number of start points for the multi homogeneous start system. 


			##Details:
					If we have all valid partitions of the degree matrix, we can find the total number of start points. This is done by 
					multiply all entries of a partition and adding that to all other products from all other partitions. 
		*/
		unsigned long long MHomogeneous::NumStartPoints() const
		{
			unsigned long long num_start_points = 0;
			for(int ii = 0; ii < valid_partitions_.size(); ii++)
			{ 
    			num_start_points += NumStartPointsForPartition(valid_partitions_[ii]);
			}
			return num_start_points;
		}

		/**
			\brief Function to create the degree matrix for a multi homogeneous start system.

			## Input: 
					target_system: System that we wish to solve. Using this we can decipher what our degree matrix is by looking 
					at each function's degree corresponding to each variable group declared. 


			## Output:
					None: This is purely used inside of a constructor. 


			##Details:
					We go through all variable groups and compute the corresponding degree vector for that variable group. Appending each 
					vector of degrees to a matrix to construct the overall degree matrix. 
		*/
		void MHomogeneous::CreateDegreeMatrix(System const& target_system)
		{
			degree_matrix_ = Mat<int>::Zero(target_system.NumNaturalFunctions(),target_system.NumTotalVariableGroups());

			var_groups_ = target_system.HomVariableGroups();
			auto affine_var_groups = target_system.VariableGroups();
			//This concatenates the affine variable groups to the hom variable groups. 
			var_groups_.insert(var_groups_.end(), affine_var_groups.begin(), affine_var_groups.end());

			int col_count = 0;
			int outer_loop = 0;
			size_t var_count = 0;
			std::vector<int> zero_column_vector(target_system.Degrees(*(var_groups_.begin())).size(), 0);

			for(std::vector<VariableGroup>::iterator it = var_groups_.begin(); it != var_groups_.end(); ++it)
			{
				outer_loop++;
  				std::vector<int> degs = target_system.Degrees(*it);
				
				std::vector<size_t> temp_v;
				for(int ii = 0; ii < (*it).size(); ++ii)
				{
					temp_v.push_back(var_count);
					var_count++;
				}
				variable_cols_.push_back(temp_v);
				
  				if(degs == zero_column_vector)
  				{
  					//check for zero column in degree matrix.
  					throw std::runtime_error("zero column in degree matrix for m-homogeneous start system!");
  				}

  				for(int ii = 0; ii <= degs.size() - 1; ++ii)
  				{
  					degree_matrix_(ii,col_count) = degs[ii];
  				}
  				col_count++;
			}

			//check for zero row in degree matrix. 
			Vec<int> zero_row_vector = degree_matrix_.row(0)*0;

			for(int ii = 0; ii < degree_matrix_.rows(); ii++)
			{

				if(degree_matrix_.row(ii) == zero_row_vector.transpose())//transpose to make comparison work correctly.
				{
					throw std::runtime_error("zero row in degree matrix for m-homogeneous start system!");
				}
			}




		}

		/**
			\brief Function to find all valid partitions inside of a degree matrix. 

			## Input: 
				target_system: System that we wish to solve. Using this we can know total number of variables, functions, and variable groups,
				 by using the systems member functions. 


			## Output:
					None: This is purely used inside of a constructor. 


			##Details:
					TODO: Fill this.
		*/
		void MHomogeneous::GenerateValidPartitions(System const& target_system)
		{
			int row = 0;
			int old_current_part_row = -1;
			int bad_choice = 0;
			Vec<int> current_partition = -1*Vec<int>::Ones(target_system.NumNaturalFunctions());
			Vec<int> variable_group_counter = Vec<int>::Zero(target_system.NumTotalVariableGroups());

			auto size_of_each_var_gp = target_system.VariableGroupSizes(); //K
			




			for(int ii = 0; ii < target_system.NumTotalVariableGroups(); ++ii)
			{
				variable_group_counter[ii] = size_of_each_var_gp[ii];
			}			    
			// std::cout << "variable_group_counter is " << std::endl;
			// std::cout << variable_group_counter << std::endl;		
			  while (row > -1)  // Algorithm will move up and down rows, kicking out to row=-1 at end
			  {
//			  	 std::cout << "current_partition before is " << std::endl;
//			  	 std::cout <<  current_partition << std::endl;
			    old_current_part_row = current_partition[row];  //Hang on to previous choice of column for this row, in case we are done with this row.
			    current_partition[row] = ChooseColumnInRow(target_system,variable_group_counter,row,current_partition[row]);  //Pick next column (var gp) for the current row (func)
//				  std::cout << "current_partition after is " << std::endl;
//				  std::cout <<  current_partition << std::endl;

			    //ChooseColumnInRow() will make this happen if it runs into col being equal to system.NumVariables()!
			    if (current_partition[row] == target_system.NumTotalVariableGroups()) // means we have exhausted all good columns for the current row, so we go back up a row
			    { //no choices for current row
			      row = row - 1;  //go back up a row
			      bad_choice = 1;
			    }
			    else  //found a good choice of column for this row!
			    {
			      	row = row + 1;  //move on to next row!
			      	if (row < target_system.NumNaturalFunctions())
			        	current_partition[row] = -1; //This allows us to consider all possible columns from left to right.
			        	//since we are starting a new row, we start with the left-most entry (ChooseColumnInRow() first increments col)
			    }
			     
			    if((row == target_system.NumNaturalFunctions()) && (!bad_choice))
			    {
			    	// std::cout << "Good partition!!!!" << std::endl;
			    	// std::cout << current_partition << std::endl;
			    	valid_partitions_.push_back(current_partition);
			     	row = row - 1; //put the counter back on the last row to try to move to the next column
			    }
			    bad_choice=0;
			  }
		}


		/**
			\brief Function to find a valid column for a given row. 

			## Input: 
				target_system: System that we wish to solve. Using this we can know total number of variables, functions, and variable groups,
				 by using the systems member functions. 
				variable_group_counter: This vector of integers will keep track of how many times we can pick from a column.
				row: This is the row in GenerateValidPartitions(). 
				column: This comes from current_partition[row] in GenerateValidPartitions(). Because of this we will walk through the 
				degree matrix from left to right.  


			## Output:
					int col: Returns the column that we have chosen for a given row. 


			##Details:
					Given a row, we search for a valid column from left to right. Since column is defaulted to -1 we will start with 
					col = 0. Variaable_group_counter will be incremented or decremented if we are coming off a good partition, or have found a 
					good column. 

		*/
		int MHomogeneous::ChooseColumnInRow(System const& target_system,Vec<int>& variable_group_counter, int row, int column)
		{
			int col = column + 1;  //We assume the current column is done and we need to increment by at least one. 
			int done = 0;		   //Note: we started at -1 so this is ok

			if (col - 1 > -1)
			{ //If we are coming off of a good partition, we need to remember to increment var_gp_ctr, 
			  //this holds how many we have chosen in a column as we move away from that column.
			  variable_group_counter[col - 1] = variable_group_counter[col - 1] + 1;// var_gp_ctr is incremented because if(var_gp_ctr[col] == 0) -> bad choice
			} 
			while (!done)
			{
				/*We have reached the end of the degree matrix. 
				Return and (current_partition[row] == target_system.NumTotalVariableGroups()) in GenerateValidPartitions() gets executed.
			    */
			    if (col == target_system.NumTotalVariableGroups())
			    {
			    	done = 1; //got to the right end of the degree matrix!
			    } 
			   
			    else
			    {
			    	/*
			    		If degree_matrix(row,col) == 0, we know that this is not a valid partition choice. We cannot pick 0 linears. 
						If variable_group_counter[col] == 0, we have picked the maximum number of choices for this column.
			    	*/
			        if ((degree_matrix_(row,col) == 0) || (variable_group_counter[col] == 0))  //bad choice, either way!
			        {	
			        	col = col + 1;
			        }
			        else  
			        {	/*means degree_matrix[row][col] > 0, so we have a possible valid choice,
			        	 col <= NumTotalVariableGroups, we have not ran out of the matrix,
			        	  and variable_group_counter[col] > 0 --> we still have choices for this column! Good column!
			            */
			            done = 1;
			            variable_group_counter[col] = variable_group_counter[col] - 1;
			        }
			    }
			}
			return col;

		}


		
		template<typename T>
		void MHomogeneous::GenerateStartPointT(Vec<T>& start_point, unsigned long long index) const
		{
			if(valid_partitions_.size() <= 0)
				throw std::runtime_error("Trying to generate MHom start points before determining valid partitions.");
			
			
			
			// First, determine which partition we are looking through
			int counter = 0;
			int partition_ii = -1;
			for (int ii = 0; ii < valid_partitions_.size(); ++ii)
			{
				counter += NumStartPointsForPartition(valid_partitions_[ii]);
				
				if(index < counter)
				{
					partition_ii = ii;
					counter -= NumStartPointsForPartition(valid_partitions_[ii]);
					index -= counter;
					break;
				}
			}
			
			if(partition_ii == -1)
				throw std::runtime_error("attempted to generate mhom start point, but index not valid");
			
			
			// Using partition ii, create dimension vector.  Then find the subscript.
			auto partition = valid_partitions_[partition_ii];
			std::vector<size_t> dim_vector(partition.size());
			for (int ii = 0; ii < partition.size(); ++ii)
			{
				dim_vector[ii] = degree_matrix_(ii, partition(ii));
			}
			
			std::vector<size_t> subscript = IndexToSubscript<size_t>(index, dim_vector);
			
			
			
			
			
			
			// Create a linear system to solve.
			size_t num_grouped_variables = NumNaturalVariables() - NumUngroupedVariables();
			Mat<T> A(partition.size(), num_grouped_variables);
			Vec<T> v(num_grouped_variables);
			Vec<T> b(partition.size());
			
			for(int ii = 0; ii < partition.size(); ++ii)
			{
				v.setZero();
				std::vector<size_t> cols = variable_cols_[partition[ii]];
				auto coeff = linprod_matrix_(ii,partition[ii])->GetCoeffs<T>(subscript[ii]);
				for(int jj = 0; jj < cols.size(); ++jj)
				{
					v(cols[jj]) = coeff[jj];
				}
				
				A.row(ii) = v;
				b(ii) = -coeff[cols.size()];
			}
			
			start_point = A.partialPivLu().solve(b);
			
			
			var_groups_[0][0]->set_current_value(start_point(0));
			var_groups_[0][1]->set_current_value(start_point(1));
			
			std::shared_ptr<node::Node> f = Function(0);
			// !!! superfluous line?
			
			

		}
		
		
		Vec<dbl> MHomogeneous::GenerateStartPoint(dbl,unsigned long long index) const
		{
			Vec<dbl> start_point(NumVariables());
			GenerateStartPointT(start_point, index);
			
			return start_point;
		}


		Vec<mpfr_complex> MHomogeneous::GenerateStartPoint(mpfr_complex,unsigned long long index) const
		{
			Vec<mpfr_complex> start_point(NumVariables());
			GenerateStartPointT(start_point, index);

			return start_point;	
		}

		inline
		MHomogeneous operator*(MHomogeneous td, std::shared_ptr<node::Node> const& n)
		{
			td *= n;
			return td;
		}


		unsigned long long MHomogeneous::NumStartPointsForPartition(Vec<int> partition) const
		{
			unsigned long long num_points = 1;
    			for(int ii = 0; ii < partition.size() ; ii++)
    			{
    				num_points *= degree_matrix_(ii,partition[ii]); 
    			}

			return num_points;
		}

	} // namespace start_system

} //namespace bertini
