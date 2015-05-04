//
//  main.cpp
//  QiPlay
//
//  Created by Collins, James B. on 3/14/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//


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
        Node(double num) {op_ = "const = " + std::to_string(num); isVar_ = true;};
        Node(int varint) {op_ = "var #" + std::to_string(varint); isVar_ = true;};
        
        std::string op_;
        std::vector<Node*> leaves;
        bool isVar_ = false;
        
        
        void addLeaf(Node* leaf);
//        void setop(std::string input) {op_ = input;}
        
        void setop(int input)
        {
            op_ = "exp^" + std::to_string(input);
        }
        
        void print()
        {
            if(!isVar_)
            {
                for(int ii = 1; ii < leaves.size(); ++ii)
                {
                    for(int jj = 0; jj < tabcount; ++jj)
                    {
                        std::cout << "\t";
                    }
                    std::cout << op_ << " leaf " << ii << ": \n";
                    tabcount++;
                    leaves[ii]->print();
                    tabcount--;
                }
            }
            else
            {
                for(int jj = 0; jj < tabcount; ++jj)
                {
                    std::cout << "\t";
                }
                std::cout << op_ << "\n";
//                tabcount--;
            }
        };
        
        void debug(std::string s)
        {
            std::cout << s;
        }
        
    };
    void Node::addLeaf(Node* leaf)
    {
        leaves.push_back(leaf);
    }
    
    
    
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
            
            sumexpr = eps[_val = phx::new_<Node>("sum")] >>
                eps[phx::bind(&Node::addLeaf,_val,new Node("null"))] >>
                multexpr[phx::bind(&Node::addLeaf,_val,_1)] >>
                *( ('-' >> multexpr[phx::bind(&Node::addLeaf,_val,_1)]) |
                  ('+' >> multexpr[phx::bind(&Node::addLeaf,_val,_1)]) );
            
            multexpr = eps[_val = phx::new_<Node>("mult")] >>
                eps[phx::bind(&Node::addLeaf,_val,new Node("null"))] >>
                subexpr[phx::bind(&Node::addLeaf,_val,_1)] >>
            *( '*' >> subexpr )[phx::bind(&Node::addLeaf,_val,_1)];
            
            subexpr = base[_val = _1] | parenexpr[_val = _1];
            
            parenexpr = eps[_val = phx::new_<Node>("exp^1")] >>
            '(' >> sumexpr[phx::bind(&Node::addLeaf,_val,new Node("null")), phx::bind(&Node::addLeaf,_val,_1)] >> ')' >>
            -( '^' >> qi::int_[phx::bind(&Node::setop,_val,_1)] );
            
            
            base = ( qi::double_[_val = phx::new_<Node>(_1)] ) |
            ( eps[_val = phx::new_<Node>("exp^1")] >>
                var[phx::bind(&Node::addLeaf,_val,new Node("null")), phx::bind(&Node::addLeaf,_val,phx::new_<Node>(_1))] >>
                -('^' >> qi::int_[phx::bind(&Node::setop,_val,_1)]) );
            
        }
        
        qi::rule<Iterator,Node*(),ascii::space_type> sumexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> multexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> subexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> parenexpr;
        qi::rule<Iterator,Node*(),ascii::space_type> base;
    };
    
    

}



int main()
{
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "\t\tPolynomial Parser\n\n";
    std::cout << "/////////////////////////////////////////////////////////\n\n";
    std::cout << "Enter your symbols separated by commas ...or [q or Q] to quit\n\n";
    
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
    client::tabcount = 0;
    while (std::getline(std::cin, str))
    {
        if (str.empty() || str[0] == 'q' || str[0] == 'Q')
            break;
        
        std::string::const_iterator iter = str.begin();
        std::string::const_iterator end = str.end();
        //[tutorial_roman_grammar_parse
        
        client::tabcount = 0;
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


