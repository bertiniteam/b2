//
//  function_tree_test.cpp
//  b2Test
//
//  Created by Collins, James B. on 4/30/15.
//  Copyright (c) 2015 West Texas A&M University. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <vector>


#include "function_tree.hpp"


#include <boost/spirit/include/qi.hpp>



#include <boost/test/unit_test.hpp>




BOOST_AUTO_TEST_SUITE(function_tree_class)

BOOST_AUTO_TEST_CASE(manual_construction_x_squared){
	
	std::shared_ptr<Variable> x = std::make_shared<Variable>("x");
	
	std::shared_ptr<Node> N = x;

	N *= N;
	
	x->set_current_value(3.5);
	

	BOOST_CHECK_EQUAL(N->Eval<dbl>() , 12.25 );
	
}


BOOST_AUTO_TEST_CASE(manual_construction_sqrt_x){
	
	std::shared_ptr<Variable> x = std::make_shared<Variable>("x");
	
	std::shared_ptr<Node> N = pow(x, 1.0/2);
	
	x->set_current_value(3.0);
	
	BOOST_CHECK_EQUAL(N->Eval<dbl>() , 1.73205080756887729352744634151);  // root 3 to 30 digits
	
}




BOOST_AUTO_TEST_SUITE_END()
