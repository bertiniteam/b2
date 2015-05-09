// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// sum_operator.h:  Declares the class SumOperator.

#ifndef b2Test_SumOperator_h
#define b2Test_SumOperator_h

#include <vector>
#include <string>   

#include "node.h"
#include "nary_operator.h"


// Node -> NaryOperator -> SumOperator
// Description: This class represents summation and difference operators.  All children are terms and are stored
// in a single vector, and a vector of bools is used to determine the sign of each term.  FreshEval method
// is defined for summation and difference.
class SumOperator : public NaryOperator
{
public:
    // Print the data for this node, as well as all it's children
    // TODO(JBC): Implement this method
    std::string PrintNode() override {return "";}
    
    
    
    //Special Behaviour: by default all terms added are positive
    virtual void AddChild(std::shared_ptr<Node> child) override
    {
        NaryOperator::AddChild(std::move(child));
        children_sign_.push_back(true);
    }
    
    
    
    //Special Behaviour: Pass bool to set sign of term: true = add, false = subtract
    void AddChild(std::shared_ptr<Node> child, bool sign)
    {
        NaryOperator::AddChild(std::move(child));
        children_sign_.push_back(sign);
    }




protected:
    // Specific implementation of FreshEval for add and subtract.
    //  If child_sign_ = true, then add, else subtract
    dbl FreshEval(dbl) override
    {
        dbl retval{};
        for(int ii = 0; ii < children_.size(); ++ii)
        {
            if(children_sign_[ii])
            {
                retval += children_[ii]->Eval<dbl>();
            }
            else
            {
                retval -= children_[ii]->Eval<dbl>();
            }
        }
        
        return retval;
    }
    
    
    
    
    mpfr FreshEval(mpfr) override
    {
        mpfr retval{};
        for(int ii = 0; ii < children_.size(); ++ii)
        {
            if(children_sign_[ii])
            {
                retval += children_[ii]->Eval<mpfr>();
            }
            else
            {
                retval -= children_[ii]->Eval<mpfr>();
            }
        }
        
        return retval;
    }

    
    
    
    
    
private:
    // Stores the sign of the particular term.  There is a one-one
    // correspondence between elements of children_sign_ and children_.  This
    // is enforced by the AddChild method below, redefined in SumOperator.
    
    // TODO(JBC): If we add method to delete child, must also delete children_sign_ entry.
    std::vector<bool> children_sign_;
    

    
    
    
    
    
    
    
    

    

};


#endif
