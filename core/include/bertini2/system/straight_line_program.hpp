//This file is part of Bertini 2.
//
//straight_line_program.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//straight_line_program.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with straight_line_program.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2021 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin eau claire

/**
\file straight_line_program.hpp 

\brief Provides the bertini::StraightLineProgram class.
*/

#ifndef BERTINI_SLP_HPP
#define BERTINI_SLP_HPP

#include <assert.h>
#include <vector>


#include "bertini2/mpfr_complex.hpp"
#include "bertini2/mpfr_extensions.hpp"
#include "bertini2/eigen_extensions.hpp"
#include "bertini2/function_tree/forward_declares.hpp"
#include "bertini2/detail/visitor.hpp"

// code copied from Bertini1's file include/bertini.h


/*

typedef struct
{
  int num_funcs;
  int num_hom_var_gp;
  int num_var_gp;
  int *type; // 0 - hom_var_gp, 1 - var_gp
  int *size; // size of the group of the user listed variables (total size = size + type)
} preproc_data;


typedef struct
{
  point_d funcVals;
  point_d parVals;
  vec_d parDer;
  mat_d Jv;
  mat_d Jp;
} eval_struct_d;

typedef struct
{
  point_mp funcVals;
  point_mp parVals;
  vec_mp parDer;
  mat_mp Jv;
  mat_mp Jp;
} eval_struct_mp;


The straight-line program structure.  This is the way that polynomials are stored internally.
typedef struct {
  int *prog;     //  The program instructions. (a big integer array)  
  int  size;     //  size of the instruction program.      
  int  memSize;  // Amount of memory it needs in workspace (for temp and final results).   
  num_t *nums;   // The array of real numbers. 
  int precision; // The precision at which evaluation should occur 
  
  // INFO NEEDED FOR M-HOM: 
  int num_var_gps;  // The total number of variable groups (i.e., m from m-hom).
  int *var_gp_sizes;  // The size of each of the groups. 
  int index_of_first_number_for_proj_trans;  // The address of the first number used in the projective transformation polynomials. 
                                                                                    
  // STOP LOCATIONS: 
  int  numInstAtEndUpdate; // instruction number at end of update. i.e. i = 0; while (i < numInstAtEndUpdate) .. 
  int  numInstAtEndParams; // instruction number at end of params. i.e. i = numInstAtEndUpdate; while (i < numInstAtEndParams) .. 
  int  numInstAtEndFnEval; // instruction number at end of function eval. i.e. i = numInstAtEndParams; while (i < numInstAtEndFnEval) .. 
  int  numInstAtEndPDeriv; // instruction number at end of param diff. i.e. i = numInstAtEndFnEval; while (i < numInstAtEndPDeriv) .. 
  int  numInstAtEndJvEval; // instruction number at end of Jv eval. i.e. i = numInstAtEndPDeriv; while (i < numInstAtEndJvEval) .. 
                           // for Jp eval: i = numInstAtEndJvEval; while (i < size) ..  

  // INPUT AMOUNTS: 
  int  numVars;  //  Number of variables in the function being computed.   
  int  numPathVars;  //  Number of path variables.  Ought to be 1 usually.   
  int  numNums;  //  Number of real numbers used in evaluation.  
  int  numConsts;  //  Number of constants.     
                                                                                   
  // OUTPUT AMOUNTS: 
  int  numPars;  //  Number of parameters   
  int  numFuncs; //  Number of coordinate functions in the homotopy.   
  int  numSubfuncs;  //  Number of subfunctions.  
                                                                                      
  // INPUT LOCATIONS: 
  int  inpVars;  //  Where the input variable values are stored.   
  int  inpPathVars;  //  Where the values of the path variables are stored.   
  int  IAddr;  //  Where the constant I is stored.  
  int  numAddr;  //  Where the first num_t is stored.  
  int  constAddr;  //  Where the first constant is stored.  
                                                                                      
  // OUTPUT LOCATIONS: 
  int  evalPars;  //  Where U(t), for given t, is stored.   
  int  evalDPars;  //  Where the derivatives of the parameters are stored.   
  int  evalFuncs;  //  Where H(x,t) is stored.   
  int  evalJVars;  //  Where the Jacobian w.r.t. vars is stored.   
  int  evalJPars;  //  Where the Jacobian w.r.t. pars is stored.   
  int  evalSubs;  //  Where the subfunctions are stored 
  int  evalJSubsV;  //  Where the derivatives of the subfunctions w.r.t. vars are stored.  
  int  evalJSubsP;  //  Where the derivatives of the subfunctions w.r.t. pars are stored.  
} prog_t;
*/



namespace bertini {

  class System; // a forward declaration, solving the circular inclusion problem







	class StraightLineProgram{
  private:
    using Nd = std::shared_ptr<node::Node>;

	public:

    /**
    The constructor -- how to make a SLP from a System.
    */
    StraightLineProgram(System const & sys);

    StraightLineProgram() = default;

		template<typename Derived>
		void Eval(Eigen::MatrixBase<Derived> const& variable_values) const
		{
      // 1. copy variable values into memory locations they're supposed to go in
      CopyVariableValues(variable_values);
		}


    template<typename Derived, typename ComplexT>
    void Eval(Eigen::MatrixBase<Derived> const& variable_values, ComplexT const& time) const
    {
      // 1. copy variable values into memory locations they're supposed to go in
      CopyVariableValues(variable_values);
      CopyPathVariable(time);
    }

    template<typename T>
    Vec<T> GetFuncVals(){}

    template<typename T>
    Mat<T> GetJacobian(){}



    inline unsigned NumTotalFunctions() const{ return num_total_functions_;}

    inline unsigned NumVariables() const{ return num_variables_;}


		/**
		\brief Get the current precision of the SLP.
		*/
		inline
		unsigned precision() const
		{
			return precision_;
		}

    /*
    change the precision of the SLP.  Downsamples from the true values.
    */
		void precision(unsigned new_precision) const;


    bool HavePathVariable() const {
      throw std::runtime_error("calling unimplemented function HavePathVariable");
      return false;
    }
	private:

    template<typename Derived>
    void CopyVariableValues(Eigen::MatrixBase<Derived> const& variable_values) const{
      throw std::runtime_error("calling unimplemented function CopyVariableValues");
    }

    template<typename ComplexT>
    void CopyPathVariable(ComplexT const& time) const{
      if (!this->HavePathVariable())
        throw std::runtime_error("calling Eval with path variable, but system doesn't have one.");
      // then actually copy the path variable into where it goes in memory
    }


    unsigned precision_ = 0;
    unsigned num_total_functions_ = 0;
    unsigned num_variables_ = 0;

    std::vector<int> instructions_;
    std::tuple< std::vector<dbl_complex>, std::vector<mpfr_complex> > memory_;
	};
	

  class SLPCompiler : public VisitorBase, public Visitor<node::Node, void>{
  private:
    using Nd = std::shared_ptr<node::Node>;
    using SLP = StraightLineProgram;

    public:

      SLP Compile(System const& sys);

      virtual void Visit(node::Function const &);
      virtual void Visit(node::SumOperator const &);
      virtual void Visit(node::Node const &);


    private:

      std::map<Nd, int> locations_encountered_symbols_;
  };



} // namespace bertini









#endif // for the ifndef include guards



