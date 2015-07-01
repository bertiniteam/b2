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
// special_number.hpp:  Declares the class SpecialNumber.



#ifndef b2Test_Differential_h
#define b2Test_Differential_h

#include <memory>
#include "function_tree/node.hpp"
#include "function_tree/symbols/symbol.hpp"
#include "function_tree/symbols/number.hpp"
// #include "function_tree/symbols/variable.hpp"



namespace bertini {
    class Variable;


    // Node -> Symbol -> SpecialNumber
    // Description: This class represents constant leaves to a function tree.  FreshEval simply returns
    // the value of the constant.
    class Differential : public virtual NamedSymbol
    {
    public:
        Differential(){};


        /**
         Input shared_ptr to a Variable.
         */
        Differential(std::shared_ptr<Variable> diff_variable, std::string var_name)
        {
            differential_variable_ = diff_variable;
            name("d" + var_name);
        }





        std::shared_ptr<Variable> GetVariable() const 
        {
            return differential_variable_;
        }


        virtual void print(std::ostream & target) const override
        {
            target << name();
        }



        /**
         Differentiates a number.  Should this return the special number Zero?
         */
        virtual std::shared_ptr<Node> Differentiate() override
        {
            return std::make_shared<Number>(0.0);
        }



        virtual ~Differential() = default;




        /**
        Compute the degree with respect to a single variable.   For differentials, the degree is 0.
        */
        virtual int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
        {
            return 0;
        }

        virtual void Homogenize(std::vector< std::shared_ptr< Variable > > const& vars, std::shared_ptr<Variable> const& homvar) override
        {
            
        }
        

    protected:
        // This should never be called for a Differential.  Only for Jacobians.
        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
        {
            if(differential_variable_ == diff_variable)
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        }

        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
        {
            if(differential_variable_ == diff_variable)
            {
                return mpfr("1.0");
            }
            else
            {
                return mpfr("0.0");
            }
        }


    private:
        std::shared_ptr<Variable> differential_variable_;
    };


} // re: namespace bertini

#endif
