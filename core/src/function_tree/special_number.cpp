#include "function_tree/symbols/special_number.hpp"



namespace bertini {


	std::shared_ptr<Node> Pi()
    {
    	return std::make_shared<SpecialNumber::Pi>();
    }

    std::shared_ptr<Node> E()
    {
    	return std::make_shared<SpecialNumber::E>();
    }

    std::shared_ptr<Node> I()
    {
    	return std::make_shared<Float>("0.0","1.0");
    }


    std::shared_ptr<Node> Two()
    {
    	return std::make_shared<Integer>(2);
    }

    std::shared_ptr<Node> One()
    {
    	return std::make_shared<Integer>(1);
    }

    std::shared_ptr<Node> Zero()
    {
    	return std::make_shared<Integer>(0);
    }

}

