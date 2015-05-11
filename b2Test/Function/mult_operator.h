// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// mult_operator.h:  Declares the class MultOperator.



#ifndef b2Test_MultOperator_h
#define b2Test_MultOperator_h

#include "node.h"
#include "nary_operator.h"



// Node -> NaryOperator -> MultOperator
// Description: This class represents multiplication operator.  All children are factors and are stored
// in a vector.  FreshEval method is defined for multiplication.
class MultOperator : public NaryOperator
{
public:
    // See node.h for description
    // TODO(JBC): Implement this
    std::string PrintNode() override {return "";}


    //Special Behaviour: by default all factors are in numerator
    virtual void AddChild(std::shared_ptr<Node> child) override
    {
        NaryOperator::AddChild(std::move(child));
        children_mult_.push_back(true);
    }
    
    
    
    //Special Behaviour: Pass bool to set sign of term: true = mult, false = divide
    void AddChild(std::shared_ptr<Node> child, bool mult)
    {
        NaryOperator::AddChild(std::move(child));
        children_mult_.push_back(mult);
    }




protected:
    // Specific implementation of FreshEval for mult and divide.
    //  If child_mult_ = true, then multiply, else divide
    dbl FreshEval(dbl) override
    {
        dbl retval{1};
        for(int ii = 0; ii < children_.size(); ++ii)
        {
            if(children_mult_[ii])
            {
                retval *= children_[ii]->Eval<dbl>();
            }
            else
            {
                retval /= children_[ii]->Eval<dbl>();
            }
        }
        
        return retval;
    }
    
    
    
    
    mpfr FreshEval(mpfr) override
    {
        mpfr retval{1};
        for(int ii = 0; ii < children_.size(); ++ii)
        {
            if(children_mult_[ii])
            {
                retval *= children_[ii]->Eval<mpfr>();
            }
            else
            {
                retval /= children_[ii]->Eval<mpfr>();
            }
        }
        
        return retval;
    }
    
    
    
    
    
    
private:
    // Stores the mult/div of a factor.  There is a one-one
    // correspondence between elements of children_sign_ and children_.  This
    // is enforced by the AddChild method, redefined in MultOperator.
    
    // TODO(JBC): If we add method to delete child, must also delete children_mult_ entry.
    std::vector<bool> children_mult_;

};




#endif
