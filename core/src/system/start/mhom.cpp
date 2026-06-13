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
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// karleigh cameron, colorado state university
// Tim Hodges, Colorado State University
// silviana amethyst, university of wisconsin eau claire

#include "bertini2/system/start/mhom.hpp"

#include "bertini2/system/blocks/block.hpp"
#include "bertini2/random.hpp"

#include <map>


BOOST_CLASS_EXPORT(bertini::start_system::MHomogeneous);


namespace bertini 
{
	using namespace bertini::node;
	
	namespace start_system 
	{

		// constructor for MHomogeneous start system, from any other *suitable* system. System can be homogeneous or non-homogeneous. 
		MHomogeneous::MHomogeneous(System const& s)
		{

			// A square multiprojective system has one equation per dimension.  Each
			// projective (homogeneous) variable group of size k spans P^{k-1}: its k
			// coordinates carry only k-1 dimensions (scale is free), so subtract one per hom
			// group.  This equals the old NumTotalFunctions()==NumVariables() for affine
			// systems and for homogenized+patched systems, but is also correct for raw systems
			// that carry projective variable groups.
			if (s.NumNaturalFunctions() != s.NumNaturalVariables() - s.NumHomVariableGroups())
				throw std::runtime_error("attempting to construct multi homogeneous start system from non-square target system");

			if (s.HavePathVariable())
				throw std::runtime_error("attempting to construct multi homogeneous start system, but target system has path varible declared already");			

			if (!s.IsPolynomial())
				throw std::runtime_error("attempting to construct multi homogeneous start system from non-polynomial target system");


			CreateDegreeMatrix(s);


			GenerateValidPartitions(s);
			
			CopyVariableStructure(s);
			
			
			// Generate the random linear-factor coefficients directly (no node::LinearProduct).
			// linear_coeffs_(ii,jj) holds degree_matrix_(ii,jj) factors over group jj's
			// variables; each factor is a row of (group_size + 1) coefficients, the trailing
			// one being the constant.  Projective groups (the leading var_groups_, indices <
			// NumHomVariableGroups) have homogeneous factors, so their constant is 0.  We
			// generate at MaxPrecisionAllowed so the block's master is precision-faithful.
			{
				auto const saved_prec = DefaultPrecision();
				DefaultPrecision(MaxPrecisionAllowed());

				linear_coeffs_ = Mat<Mat<mpfr_complex>>(degree_matrix_.rows(), degree_matrix_.cols());
				for (Eigen::Index ii = 0; ii < degree_matrix_.rows(); ++ii)
					for (Eigen::Index jj = 0; jj < degree_matrix_.cols(); ++jj)
					{
						const int d = degree_matrix_(ii, jj);
						if (d == 0)
							continue;
						const Eigen::Index gsize = static_cast<Eigen::Index>(var_groups_[jj].size());
						const bool projective = (static_cast<size_t>(jj) < s.NumHomVariableGroups());
						Mat<mpfr_complex> C(d, gsize + 1);
						for (Eigen::Index f = 0; f < d; ++f)
						{
							for (Eigen::Index k = 0; k < gsize; ++k)
								C(f, k) = mpfr_complex(mpfr_float(RandomRat()), mpfr_float(RandomRat()));
							C(f, gsize) = projective ? mpfr_complex(0)
							                         : mpfr_complex(mpfr_float(RandomRat()), mpfr_float(RandomRat()));
						}
						linear_coeffs_(ii, jj) = std::move(C);
					}

				DefaultPrecision(saved_prec);
			}

			// Homogenize the variable structure (adds a homogenizing variable per affine
			// group; projective groups already carry their own) and copy the target's patch.
			// There are no node functions to homogenize -- the start system evaluates through
			// the products-of-linears block built below.
			if (s.IsHomogeneous())
				Homogenize();

			if (s.IsPatched())
				CopyPatches(s);

			// Build the products-of-linears evaluation block from linear_coeffs_, so the start
			// system evaluates via the block (the SLP compiler cannot compile linear-product
			// node trees, which is what blocked MHom in a homotopy).  Per function we assemble
			// an augmented coefficient matrix (one row per factor, length NumVariables()+1),
			// placing each factor's coefficients by VARIABLE IDENTITY into the full variable
			// ordering.  A factor's constant multiplies the group's homogenizing variable when
			// the system is homogeneous (so it goes in that variable's column), else it is the
			// augmented constant in the trailing column.  Built at MaxPrecisionAllowed so the
			// block's master is precision-faithful.
			{
				auto const saved_prec = DefaultPrecision();
				DefaultPrecision(MaxPrecisionAllowed());
				this->precision(MaxPrecisionAllowed());

				const VariableGroup& vars = this->Variables();
				std::map<node::Node const*, Eigen::Index> col_of;
				for (Eigen::Index c = 0; c < static_cast<Eigen::Index>(vars.size()); ++c)
					col_of[vars[c].get()] = c;
				const Eigen::Index n = static_cast<Eigen::Index>(this->NumVariables());

				std::vector<Mat<mpfr_complex>> per_function;
				per_function.reserve(static_cast<size_t>(degree_matrix_.rows()));

				for (Eigen::Index ii = 0; ii < degree_matrix_.rows(); ++ii)
				{
					std::vector<Vec<mpfr_complex>> rows;
					for (Eigen::Index g = 0; g < degree_matrix_.cols(); ++g)
					{
						const int d = degree_matrix_(ii, g);
						if (d == 0)
							continue;

						const Mat<mpfr_complex>& C = linear_coeffs_(ii, g);
						const VariableGroup& gvars = var_groups_[static_cast<size_t>(g)];
						// the homogenizing variable of an affine group, if the system is
						// homogenized; projective groups (g < num_hom_groups_) have none.
						std::shared_ptr<node::Variable> hom_var;
						if (static_cast<size_t>(g) >= num_hom_groups_ &&
						    !this->HomogenizingVariables().empty())
							hom_var = this->HomogenizingVariables()[static_cast<size_t>(g) - num_hom_groups_];

						for (Eigen::Index f = 0; f < d; ++f)
						{
							Vec<mpfr_complex> row = Vec<mpfr_complex>::Zero(n + 1);
							for (size_t k = 0; k < gvars.size(); ++k)
								row(col_of.at(gvars[k].get())) = C(f, static_cast<Eigen::Index>(k));
							const mpfr_complex& constant = C(f, static_cast<Eigen::Index>(gvars.size()));
							if (hom_var && col_of.count(hom_var.get()))
								row(col_of.at(hom_var.get())) = constant;
							else
								row(n) = constant;   // augmented constant (0 for projective groups)
							rows.push_back(std::move(row));
						}
					}

					Mat<mpfr_complex> M(static_cast<Eigen::Index>(rows.size()), n + 1);
					for (size_t r = 0; r < rows.size(); ++r)
						M.row(static_cast<Eigen::Index>(r)) = rows[r].transpose();
					per_function.push_back(std::move(M));
				}

				this->AddBlock(blocks::ProductsOfLinearsBlock(static_cast<size_t>(n), std::move(per_function)));

				DefaultPrecision(saved_prec);
				this->precision(saved_prec);
			}

		}// M-Hom constructor

		
		MHomogeneous& MHomogeneous::operator*=(Nd const& n)
		{
			System::operator*=(n);
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
			for(size_t ii = 0; ii < valid_partitions_.size(); ii++)
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
			num_hom_groups_ = var_groups_.size();   // the leading var_groups_ entries are projective
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
				for(size_t ii = 0; ii < (*it).size(); ++ii)
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

  				for(size_t ii = 0; ii < degs.size(); ++ii)
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
			int bad_choice = 0;
			Vec<int> current_partition = -1*Vec<int>::Ones(target_system.NumNaturalFunctions());
			Vec<int> variable_group_counter = Vec<int>::Zero(target_system.NumTotalVariableGroups());

			// Capacity per group = the number of functions that group may be assigned, which
			// is the group's declared dimension.  Use var_groups_ (the groups as concatenated
			// for the degree matrix in CreateDegreeMatrix: hom groups then affine groups),
			// NOT target_system.VariableGroupSizes(): once the target is homogenized the
			// latter counts the added homogenizing variable, doubling the capacity and
			// admitting invalid over-filled partitions (e.g. both functions in one size-1
			// group), whose linear solve is singular and yields NaN start points.
			for(size_t ii = 0; ii < target_system.NumTotalVariableGroups(); ++ii)
			{
				// A projective (homogeneous) group of size k spans P^{k-1}, so it takes only
				// k-1 functions; an affine group of size m takes m.  The leading num_hom_groups_
				// entries of var_groups_ are the projective ones.
				const int dim = static_cast<int>(var_groups_[ii].size()) - (ii < num_hom_groups_ ? 1 : 0);
				variable_group_counter[ii] = dim;
			}
			// std::cout << "variable_group_counter is " << std::endl;
			// std::cout << variable_group_counter << std::endl;		
			  while (row > -1)  // Algorithm will move up and down rows, kicking out to row=-1 at end
			  {
//			  	 std::cout << "current_partition before is " << std::endl;
//			  	 std::cout <<  current_partition << std::endl;
			    current_partition[row] = ChooseColumnInRow(target_system,variable_group_counter,row,current_partition[row]);  //Pick next column (var gp) for the current row (func)
//				  std::cout << "current_partition after is " << std::endl;
//				  std::cout <<  current_partition << std::endl;

			    //ChooseColumnInRow() will make this happen if it runs into col being equal to system.NumVariables()!
			    if (current_partition[row] == static_cast<int>(target_system.NumTotalVariableGroups())) // means we have exhausted all good columns for the current row, so we go back up a row
			    { //no choices for current row
			      row = row - 1;  //go back up a row
			      bad_choice = 1;
			    }
			    else  //found a good choice of column for this row!
			    {
			      	row = row + 1;  //move on to next row!
			      	if (row < static_cast<int>(target_system.NumNaturalFunctions()))
			        	current_partition[row] = -1; //This allows us to consider all possible columns from left to right.
			        	//since we are starting a new row, we start with the left-most entry (ChooseColumnInRow() first increments col)
			    }
			     
			    if((row == static_cast<int>(target_system.NumNaturalFunctions())) && (!bad_choice))
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
			    if (col == static_cast<int>(target_system.NumTotalVariableGroups()))
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
			unsigned long long counter = 0;
			int partition_ii = -1;
			for (size_t ii = 0; ii < valid_partitions_.size(); ++ii)
			{
				counter += NumStartPointsForPartition(valid_partitions_[ii]);
				
				if(index < counter)
				{
					partition_ii = static_cast<int>(ii);
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
			
			
			
			
			
			
			// Create the linear system whose solution is this start point, in natural
			// coordinates.  Each function assigned by the partition contributes its chosen
			// linear factor = 0:
			//   affine group (m vars):   m factors -> m rows; the factor's trailing coeff is
			//                            its constant, so it moves to the right-hand side.
			//   projective group (k coords): k-1 factors -> k-1 rows; these are homogeneous
			//                            (trailing coeff 0), so the per-group block is k-1 by
			//                            k -- one short of square.  We pin one coordinate of
			//                            each projective group to 1 (an affine-chart normaliza-
			//                            tion) to pick a representative; HomogenizePoint then
			//                            rescales it onto the patch.  Adding one such row per
			//                            projective group makes A exactly square.
			size_t num_grouped_variables = NumNaturalVariables() - NumUngroupedVariables();
			Mat<T> A = Mat<T>::Zero(static_cast<Eigen::Index>(num_grouped_variables),
			                        static_cast<Eigen::Index>(num_grouped_variables));
			Vec<T> b = Vec<T>::Zero(static_cast<Eigen::Index>(num_grouped_variables));

			// coefficients come from linear_coeffs_ (mpfr master); cast each to the working
			// type T (a no-op widen for mpfr, a narrowing for dbl).
			auto as_T = [](mpfr_complex const& z) -> T {
				if constexpr (std::is_same<T, dbl>::value) return dbl(z);
				else return z;
			};
			for(int ii = 0; ii < partition.size(); ++ii)
			{
				std::vector<size_t> cols = variable_cols_[partition[ii]];
				const Mat<mpfr_complex>& C = linear_coeffs_(ii, partition[ii]);
				const Eigen::Index f = static_cast<Eigen::Index>(subscript[ii]);
				for(size_t jj = 0; jj < cols.size(); ++jj)
				{
					A(ii, static_cast<Eigen::Index>(cols[jj])) = as_T(C(f, static_cast<Eigen::Index>(jj)));
				}
				b(ii) = -as_T(C(f, static_cast<Eigen::Index>(cols.size())));   // affine: the constant term; projective: 0
			}

			// normalization row per projective group: pin its last coordinate to 1.
			Eigen::Index extra = static_cast<Eigen::Index>(partition.size());
			for (size_t g = 0; g < num_hom_groups_; ++g)
			{
				A(extra, static_cast<Eigen::Index>(variable_cols_[g].back())) = T(1);
				b(extra) = T(1);
				++extra;
			}

			Vec<T> affine_solution = A.partialPivLu().solve(b);

			// The linear solve gives the start point in affine (dehomogenized) coordinates;
			// lift it onto the homogenized + patched coordinate system the homotopy is
			// tracked in (HomogenizePoint inserts the homogenizing coordinates per affine
			// group and rescales onto the patch, and is a no-op if not homogenized).
			start_point = this->HomogenizePoint(affine_solution);

			// Return the point at the current working precision.  The linear-factor
			// coefficients and the patch are stored at MaxPrecisionAllowed, which would
			// otherwise leak into the start point and mismatch the tracker's working
			// precision (the adaptive tracker begins at the ambient precision).
			if constexpr (!std::is_same<T, dbl>::value)
				for (Eigen::Index i = 0; i < start_point.size(); ++i)
					start_point(i).precision(DefaultPrecision());
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
