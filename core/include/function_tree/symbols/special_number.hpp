//This file is part of Bertini 2.0.
//
//special_number.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//special_number.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with special_number.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Created by Collins, James B. on 4/30/15.
//
//
// special_number.hpp:  Declares the class SpecialNumber.



#ifndef b2Test_SpecialNumber_h
#define b2Test_SpecialNumber_h


#include "function_tree/symbols/symbol.hpp"

#include "function_tree/operators/arithmetic.hpp"
#include "function_tree/operators/trig.hpp"

#include <cmath>



namespace bertini {
    
    using ::acos;
    using ::exp;

    
    namespace SpecialNumber{

	    class Pi : public NamedSymbol
	    {
	    public:
	    	Pi() : NamedSymbol("pi")
	    	{}

	    	virtual ~Pi() = default;



	    	void print(std::ostream & target) const override
			{
				target << name();
			}


			/**
	         Differentiates a number.  Should this return the special number Zero?
	         */
	        std::shared_ptr<Node> Differentiate() override
	        {
	            return std::make_shared<Float>(0.0);
	        }

	        /**
			Compute the degree with respect to a single variable.

			For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		    */
		    int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		    {
		    	return 0;
		    }


		    int Degree(VariableGroup const& vars) const override
			{
				return 0;
			}

			std::vector<int> MultiDegree(VariableGroup const& vars) const override
			{
				return std::vector<int>(vars.size(), 0);
			}


		    void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
			{
				
			}

			bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
			{
				return true;
			}
			
			/**
			Check for homogeneity, with respect to a variable group.
			*/
			bool IsHomogeneous(VariableGroup const& vars) const override
			{
				return true;
			}


			/**
			 Change the precision of this variable-precision tree node.
			 
			 \param prec the number of digits to change precision to.
			 */
			virtual void precision(unsigned int prec) override
			{
				auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
				val_pair.first.precision(prec);
			}



			
	    private:
	    	// Return value of constant
	        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
	        {
	            return acos(-1.0);
	        }

	        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
	        {
	            return mpfr(acos(boost::multiprecision::mpfr_float("-1.0")),0.0);
	        }

	    };




	    class E : public NamedSymbol
	    {
	    public:
	    	E() : NamedSymbol("e")
	    	{}

	    	virtual ~E() = default;



	    	void print(std::ostream & target) const override
			{
				target << name();
			}


			/**
	         Differentiates a number.  Should this return the special number Zero?
	         */
	        std::shared_ptr<Node> Differentiate() override
	        {
	            return std::make_shared<Float>(0.0);
	        }

	        /**
			Compute the degree with respect to a single variable.

			For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
		    */
		    int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
		    {
		    	return 0;
		    }


		    int Degree(VariableGroup const& vars) const override
			{
				return 0;
			}

			std::vector<int> MultiDegree(VariableGroup const& vars) const override
			{
				return std::vector<int>(vars.size(), 0);
			}


		    void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
			{
				
			}

			bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
			{
				return true;
			}
			
			/**
			Check for homogeneity, with respect to a variable group.
			*/
			bool IsHomogeneous(VariableGroup const& vars) const override
			{
				return true;
			}

			/**
			 Change the precision of this variable-precision tree node.
			 
			 \param prec the number of digits to change precision to.
			 */
			virtual void precision(unsigned int prec) override
			{
				auto& val_pair = std::get< std::pair<mpfr,bool> >(current_value_);
				val_pair.first.precision(prec);
			}

	    private:
	    	// Return value of constant
	        dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
	        {
	            return dbl(exp(1.0f),0.0);
	        }

	        mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
	        {
	            return mpfr(exp(boost::multiprecision::mpfr_float("1.0")),0.0);
	        }


	    };

	} // special number namespace

	std::shared_ptr<Node> Pi();

	std::shared_ptr<Node> E();

	std::shared_ptr<Node> I();

    std::shared_ptr<Node> Two();

    std::shared_ptr<Node> One();

    std::shared_ptr<Node> Zero();



} // re: namespace bertini
    
#endif





//     // Node -> Symbol -> SpecialNumber
//     // Description: This class represents constant leaves to a function tree.  FreshEval simply returns
//     // the value of the constant.
//     class SpecialNumber : public virtual NamedSymbol
//     {
//     public:
//         SpecialNumber(){};
        
        
//         /**
//          Input string num denotes the type of special number this node is.
//          "pi" or "Pi" is the number pi (Stored to 200 digits and obtained from Maple)
//          "e" or "E" is exp(1) (Stored to 200 digits and obtained from Maple)
//          "i" or "I" or "1i" is the imaginary number
//          */
//         SpecialNumber(std::string snum)
//         {
//             using namespace boost::math::constants;
            
//             if(snum == "pi" || snum == "Pi")
//             {
                
//                 std::get< std::pair<dbl,bool> >(current_value_).first = dbl(3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303820,0.0);
//                 std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303820", "0.0");
//                 name("pi");
//             }
//             else if(snum == "e" || snum == "E")
//             {
//                 std::get< std::pair<dbl,bool> >(current_value_).first = dbl(2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190,0.0);
//                 std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr("2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190", "0.0");
//                 name("e");
//             }
//             else if(snum == "i" || snum == "I" || snum == "1i")
//             {
//                 std::get< std::pair<dbl,bool> >(current_value_).first = dbl(0.0,1.0);
//                 std::get< std::pair<mpfr,bool> >(current_value_).first = mpfr("0.0","1.0");
//                 name("i");
//             }
//         }
        
        
        

        
    
// 	    void print(std::ostream & target) const override
// 	    {
// 	    target << std::get< std::pair<mpfr,bool> >(current_value_).first;
// 	    }
	    
// 	    void Reset() override
// 	    {
// 	        // nothing to reset here
// 	    }
	        
	        
// 	    /**
// 	     Differentiates a number.  Should this return the special number Zero?
// 	     */
// 	    std::shared_ptr<Node> Differentiate() override
// 	    {
// 	        return std::make_shared<Float>(0.0);
// 	    }

    


//         /**
//         Compute the degree with respect to a single variable.

//         For transcendental functions, the degree is 0 if the argument is constant, otherwise it's undefined, and we return -1.
//         */
//         int Degree(std::shared_ptr<Variable> const& v = nullptr) const override
//         {
//             return 0;
//         }

//         int Degree(VariableGroup const& vars) const override
// 		{
// 			return 0;
// 		}

// 		std::vector<int> MultiDegree(VariableGroup const& vars) const override
// 		{
// 			return std::vector<int>(vars.size(), 0);
// 		}

//         void Homogenize(VariableGroup const& vars, std::shared_ptr<Variable> const& homvar) override
//         {
            
//         }

//         bool IsHomogeneous(std::shared_ptr<Variable> const& v = nullptr) const override
// 		{
// 			return true;
// 		}

// 		/**
// 		Check for homogeneity, with respect to a variable group.
// 		*/
// 		bool IsHomogeneous(VariableGroup const& vars) const override
// 		{
// 			return true;
// 		}
		
//     virtual ~SpecialNumber() = default;
    
    
// protected:
//     // Return value of constant
//     dbl FreshEval(dbl, std::shared_ptr<Variable> diff_variable) override
//     {
//         return std::get< std::pair<dbl,bool> >(current_value_).first;
//     }
    
//     mpfr FreshEval(mpfr, std::shared_ptr<Variable> diff_variable) override
//     {
//         return std::get< std::pair<mpfr,bool> >(current_value_).first;
//     }

    
// };


