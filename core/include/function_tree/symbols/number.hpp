//This file is part of Bertini 2.0.
//
//number.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//number.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with number.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// number.hpp:  Declares the class Number.



#ifndef b2Test_Number_h
#define b2Test_Number_h


#include "function_tree/symbols/symbol.hpp"


namespace bertini {



	// Node -> Symbol -> Number
	// Description: This class represents constant leaves to a function tree.  FreshEval simply returns
	// the value of the constant.
	class Number : public virtual Symbol
	{
	public:
		Number(){};


		Number(int val)
		{
            std::get< std::pair<dbl,bool> >(current_value_).first = dbl(val,0.0);
            std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(val,0.0);
		}

		Number(double val)
		{
            std::get< std::pair<dbl,bool> >(current_value_).first = dbl(val,0.0);
            std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(val,0.0);
		}

        Number(double rval, double ival)
        {
            std::get< std::pair<dbl,bool> >(current_value_).first = dbl(rval,ival);
            std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(rval,ival);
		}



		// Ctor that reads in a string for the real component and converts it to a number of the appropriate type
		Number(std::string sval)
		{
            std::get< std::pair<dbl,bool> >(current_value_).first = dbl(stod(sval), 0.0);
			std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(sval, "0.0");
		}

        // Ctor that reads in two strings for a complex number and converts it to a number of the appropriate type
        Number(std::string srval, std::string sival)
        {
            std::get< std::pair<dbl,bool> >(current_value_).first = dbl(stod(srval), stod(sival));
            std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr(srval, sival);
        }






		void print(std::ostream & target) const override
		{
			target << std::get< std::pair<mpfr,bool> >(current_value_).first;
		}

		void Reset() override
		{
			// nothing to reset here
		}


        /**
         Differentiates a number.  Should this return the special number Zero?
         */
        std::shared_ptr<Node> Differentiate() override
        {
            return std::make_shared<Number>(0.0);
        }

       


		/**
		Compute the degree with respect to a single variable.

		For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
	    */
	    int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
	    {
	    	return 0;
	    }


	    int Degree(std::vector< std::shared_ptr<Variable > > const& vars) const override
		{
			return 0;
		}

		std::vector<int> MultiDegree(std::vector< std::shared_ptr<Variable> > const& vars) const override
		{
			return std::vector<int>(vars.size(), 0);
		}


	    void Homogenize(std::vector< std::shared_ptr< Variable > > const& vars, std::shared_ptr<Variable> const& homvar) override
		{
			
		}

		bool IsHomogeneous() const override
		{
			return true;
		}
		
		virtual ~Number() = default;


	protected:
		// Return value of constant
        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            return std::get< std::pair<dbl,bool> >(current_value_).first;
        }

        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            return std::get< std::pair<mpfr,bool> >(current_value_).first;
        }

	};


} // re: namespace bertini

#endif
