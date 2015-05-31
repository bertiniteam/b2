
// the purpose of this file is the three non-comment lines.  it has no other purpose than to provide a place for them.



//TODO: make the DYN_LINK change depending on the targeted architecture.  some need it, others don't.
//if used, this BOOST_TEST_DYN_LINK appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_DYN_LINK

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "Bertini 2 Class Testing"
#include <boost/test/unit_test.hpp>




// the bottom of this file is intentionally blank.  this is the 'main' .cpp file for the built boost unit test suite for bertini 2 classes.
