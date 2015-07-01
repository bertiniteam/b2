//This file is part of Bertini 2.0.
//
//function.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//function.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with function.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 6/11/15.
//
//
// jacobian.hpp:  Declares the class Jacobian.


#ifndef b2Test_Jacobian_h
#define b2Test_Jacobian_h


#include "function_tree/node.hpp"
#include "function_tree/symbols/variable.hpp"



namespace bertini {
    
    
    /**
     Node -> Function
     This class defines a function.  It stores the entry node for a particular functions tree.
     */
    class Jacobian : public Function
    {
    public:
        
        
        /**
         The default constructor
         */
        Jacobian()
        {};
        
        
        /**
         Constructor defines entry node at construct time.
         */
        Jacobian(const std::shared_ptr<Node> & entry) : Function(entry)
        {
        }
        
        
        /**
         Jacobians must be evaluated with EvalJ, so that when current_diff_variable changes
         the Jacobian is reevaluated.
         */
        template<typename T>
        T Eval(std::shared_ptr<Variable> diff_variable = nullptr) = delete;
        
        
        // Evaluate the node.  If flag false, just return value, if flag true
        //  run the specific FreshEval of the node, then set flag to false.
        //
        // Template type is type of value you want returned.
        template<typename T>
        T EvalJ(std::shared_ptr<Variable> diff_variable)
        {
            //		this->print(std::cout);
            //		std::cout << " has value ";
            
            auto& val_pair = std::get< std::pair<T,bool> >(current_value_);
            auto& variable_pair = std::get< std::pair<T,std::shared_ptr<Variable>> >(current_diff_variable_);
            if(diff_variable != variable_pair.second)
            {
                //            std::cout << "Fresh Eval\n";
                variable_pair.second = diff_variable;
                Reset();
                T input{};
                val_pair.first = FreshEval(input, diff_variable);
                val_pair.second = true;
            }
            
            
            //		std::cout << val_pair.first << std::endl;
            return val_pair.first;
        }
        



        

        
        virtual ~Jacobian() = default;
        
      
        
        

        
        
    private:
        std::tuple< std::pair<dbl,std::shared_ptr<Variable>>, std::pair<mpfr,std::shared_ptr<Variable>> > current_diff_variable_;
        
        /**
        The computation of degree for Jacobians is challenging and not correctly implemented, so it is private.
        */
        int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
        {
            return entry_node_->Degree(v);
        }

        /**
        The computation of degree for Jacobians is challenging and not correctly implemented, so it is private.
        */
        int Degree(std::vector< std::shared_ptr<Variable > > const& vars) const override
		{
			return entry_node_->Degree(vars);
		}

		/**
		 Compute the multidegree with respect to a variable group.  This is for homogenization, and testing for homogeneity.  
	    */
		std::vector<int> MultiDegree(std::vector< std::shared_ptr<Variable> > const& vars) const override
		{
			
			std::vector<int> deg(vars.size());
			for (auto iter = vars.begin(); iter!= vars.end(); ++iter)
			{
				*(deg.begin()+(iter-vars.begin())) = this->Degree(*iter);
			}
			return deg;
		}

    };
    
    
} // re: namespace bertini



#endif
