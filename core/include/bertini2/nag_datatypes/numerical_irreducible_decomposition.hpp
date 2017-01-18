//This file is part of Bertini 2.
//
//bertini2/nag_datatypes/numerical_irreducible_decomposition.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_datatypes/numerical_irreducible_decomposition.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_datatypes/numerical_irreducible_decomposition.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/nag_datatypes/numerical_irreducible_decomposition.hpp 

\brief Provides the NumericalIrreducibleDecomposition data type for Bertini2.
*/


#pragma once

#include "bertini2/nag_datatypes/witness_set.hpp"

namespace bertini {

	namespace nag_datatype {


		/**
		\brief The NID data type for Bertini2.

		*/
		template<typename ComplexT, typename SystemT = System, template<typename> class ObjManagementP = policy::Copy >
		class NumericalIrreducibleDecomposition 
		{
			using WS = WitnessSet<ComplexT, SystemT, ObjManagementP>;
			using WSCont = std::vector<WS>;

			WSCont finished_witness_sets_;

		public:
			std::vector<int> NonEmptyCodimensions() const
			{}

			const WSCont& GetWitnessSets() const
			{
				return finished_witness_sets_;
			}

			WSCont WitnessSetsOfDim(int dim)
			{
				WSCont w_correct_dim;
				for (const auto& w : finished_witness_sets_)
					if (w.Dimension() == dim)
						w_correct_dim.push_back(w);

				return w_correct_dim;
			} 

		};


	} // nag_datatype

}//namespace bertini