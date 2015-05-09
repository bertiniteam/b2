// BoilerPlate licence
//  Created by Collins, James B. on 4/30/15.
//
//
// input_parser.h:  Declares the class InputParser.

#ifndef b2Test_InputParser_h
#define b2Test_InputParser_h

#ifndef b2Test_Grammar_h
#define b2Test_Grammar_h

#define BOOST_RESULT_OF_USE_DECLTYPE
#define BOOST_SPIRIT_USE_PHOENIX_V3 1

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/bind.hpp>



#include "Node.h"
#include "sum_operator.h"
#include "mult_operator.h"
#include "negate_operator.h"
#include "exp_operator.h"
#include "constant.h"
#include "variable.h"





// InputParser
// Description: This class will hold all grammars and methods needed to parse an
// input file for Bertini.
class InputParser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
public:
    
private:
    //
    void addVars(std::vector<char> const& c)
    {
        std::string str;
        for(auto it : c)
        {
            str = str + it;
        }
        var.add(str,varcount);
        varcount++;
    }

    
    
    
    // This struct defines a new parsing type for qi which
    //  contains the variable names input from the file
    struct var_ : qi::symbols<char,int> var;
    
    // This is used when parsing the variable list in the
    //  input file.
    int varcount = 0;
};





#endif
