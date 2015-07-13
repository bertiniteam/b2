#include "limbo.hpp"



namespace bertini{

  std::vector<size_t> IndexToSubscript(size_t index, std::vector<size_t> dimensions)
  {

    std::vector< size_t > subscripts(dimensions.size());//for forming a subscript from an index

    std::vector<size_t> k(dimensions.size(),1);
    for (int ii = 0; ii < dimensions.size()-1; ++ii)
      k[ii+1] = k[ii]*dimensions[ii];


    if (index >= k.back()*dimensions.back())
      throw std::out_of_range("in IndexToSubscript, index exceeds max based on dimension sizes");


    for (int ii = dimensions.size()-1; ii >= 0; --ii)
    {
      size_t I = index%k[ii];
      size_t J = (index - I) / k[ii];
      subscripts[ii] = J;
      index = I;
    }

    return subscripts;
  }


} // re: namespace bertini

