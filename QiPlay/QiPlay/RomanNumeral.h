//
//  RomanNumeral.h
//  QiPlay
//
//  Created by Collins, James B. on 3/15/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#ifndef QiPlay_RomanNumeral_h
#define QiPlay_RomanNumeral_h
#include <boost/spirit/include/qi.hpp>


using namespace boost::spirit;

struct hundreds_ : qi::symbols<char,int>
{
    hundreds_()
    {
        add("C",100);
        add("CC",200);
        add("CCC",300);
        add("CD",400);
        add("D",500);
        add("DC",600);
        add("DCC",700);
        add("DCCC",800);
        add("CM",900);
    }
} hundreds_;

struct tens_ : qi::symbols<char,int>
{
    tens_()
    {
        add("X",10);
        add("XX",20);
        add("XXX",30);
        add("XL",40);
        add("L",50);
        add("LX",60);
        add("LXX",70);
        add("LXXX",80);
        add("XC",90);
    }
} tens_;


struct ones_ : qi::symbols<char,int>
{
    ones_()
    {
        add("I",1);
        add("II",2);
        add("III",3);
        add("IV",4);
        add("V",5);
        add("VI",6);
        add("VII",7);
        add("VIII",8);
        add("IX",9);
    }
} ones_;



#endif
