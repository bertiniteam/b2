//
//  main.cpp
//  QiPlay
//
//  Created by Collins, James B. on 3/14/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//
#define BOOST_SPIRIT_USE_PHOENIX_V3 1


#include <string>
#include <iostream>
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix/object/new.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>






namespace client
{
    
    
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    struct var_ : qi::symbols<char,int>
    {
        
    }var;
    
    void addPlus()
    {
        std::cout << "Create a plus operator\n";
    }
    
    void addVar(int varint)
    {
        std::cout << "Create a variable object " << varint << std::endl;
    }
    
    
    int varcount = 0;
    int tabcount = 0;
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

    template<typename Iterator>
    struct variables : qi::grammar<Iterator,ascii::space_type>
    {
        
        variables() : variables::base_type(start)
        {
            namespace phx = boost::phoenix;
            using qi::_1;
            
            start = (*qi::alnum)[&addVars] % ",";
        }
        
        qi::rule<Iterator,ascii::space_type> start;
    };
    
    
    class Node
    {
    public:
        Node() {op_ = "not set";};
        Node(std::string op) {op_ = op; isVar_ = false;};
        Node(std::string op, bool isVar) {op_ = op; isVar_ = isVar;};
        
        std::string op_;
        std::vector<Node*> leaves;
        bool isVar_ = false;
        
        void setop(std::string input) {op_ = input;};
        void print()
        {
            if(!isVar_)
            {
                for(int ii = 1; ii < leaves.size(); ++ii)
                {
                    std::cout << op_ << "leaf " << ii << ": \n";
//                    tabcount++;
//                    for(int jj = 0; jj < tabcount; ++jj)
//                    {
//                        std::cout << "\t";
//                    }
                    leaves[ii]->print();
//                    tabcount--;
                }
            }
            else
            {
                std::cout << op_ << "\n";
                tabcount--;
            }
        };
        void addLeaf(Node* leaf)
        {
            leaves.push_back(leaf);
        }
        
    };
    
    
    
    template<typename Iterator>
    struct polynomial : qi::grammar<Iterator,Node*(),ascii::space_type>
    {
        
        polynomial() : polynomial::base_type(sumexpr)
        {
            namespace phx = boost::phoenix;
            using qi::_1;
            using qi::_val;
            using qi::eps;
            
            
//            sumexpr = eps[_val = phx::new_<Node>("add")] >>
//            eps[phx::bind(&Node::addLeaf,_val,new Node("null",true))] >>
//            var[phx::bind(&Node::print,_val)] >>
//            *( ('-' >> var) |
//              ('+' >> var) );
            
            sumexpr = eps[_val = phx::new_<Node>("add")] >>
                eps[phx::bind(&Node::addLeaf,_val,new Node("null"))] >>
                multexpr[phx::bind(&Node::addLeaf,_val,_1)] >>
                *( ('-' >> multexpr[phx::bind(&Node::addLeaf,_val,_1)]) |
                  ('+' >> multexpr[phx::bind(&Node::addLeaf,_val,_1)]) );
            
            multexpr = eps[_val = phx::new_<Node>("mult"),phx::bind(&Node::print,_val)] >>
                eps[phx::bind(&Node::addLeaf,_val,new Node("null"))] >>
                subexpr[phx::bind(&Node::addLeaf,_val,_1)] >> eps[phx::bind(&Node::print,_val)] >>
                *( '*' >> multexpr )[phx::bind(&Node::addLeaf,_val,_1)];
            
            subexpr = monomial[_val = phx::new_<Node>("varmon",true)] | parenexpr[_val = phx::new_<Node>("varpar",true)];
            
            parenexpr = '(' >> sumexpr >> ')' >> -('^' >> qi::int_);
            
            monomial = base >> *('*' >> base);
            
            base = qi::double_ | (var >> -('^' >> qi::int_));
            
        }
        
        qi::rule<Iterator,Node*(),ascii::space_type> sumexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> multexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> subexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> parenexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> monomial;
        qi::rule<Iterator,ascii::space_type> base;
    };
    
    

}



int main()
{
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "\t\tPolynomial Parser\n\n";
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Enter your symbols ...or [q or Q] to quit\n\n";
    
    typedef std::string::const_iterator iterator_type;
    typedef client::variables<iterator_type> variables;
    typedef client::polynomial<iterator_type> polys;
    
    variables var_parser; // Our grammar
    
    std::string str;
    std::getline(std::cin, str);
    iterator_type it = str.begin();
    iterator_type end = str.end();
    phrase_parse(it,end,var_parser,client::ascii::space);
    
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Enter your polynomial ...or [q or Q] to quit\n\n";
    
    polys poly_parser;
    client::Node* test;
    while (std::getline(std::cin, str))
    {
        if (str.empty() || str[0] == 'q' || str[0] == 'Q')
            break;
        
        std::string::const_iterator iter = str.begin();
        std::string::const_iterator end = str.end();
        //[tutorial_roman_grammar_parse
        
        test = new client::Node();
        bool r = phrase_parse(iter, end, poly_parser,client::ascii::space,test);
        if(test != 0)
        {
            test->print();
        }
        else{
            std::cout << "didn't work\n";
        }
        
        
        if (r && iter == end)
        {
            std::cout << "-------------------------\n";
            std::cout << "Parsing succeeded\n";
//            std::cout << "result = " << result << std::endl;
            std::cout << "-------------------------\n";
        }
        else
        {
            std::string rest(iter, end);
            std::cout << "-------------------------\n";
            std::cout << "Parsing failed\n";
            std::cout << "stopped at: \": " << rest << "\"\n";
            std::cout << "-------------------------\n";
        }
        //]
    }
    
    std::cout << "Bye... :-) \n\n";
    return 0;
}


