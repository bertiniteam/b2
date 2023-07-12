#include "random_export.hpp"


namespace bertini{
namespace python{



void ExportRandom(){

	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".random");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("random") = new_submodule;
	scope new_submodule_scope = new_submodule;



	def("complex_in_minus_one_to_one", bertini::multiprecision::rand,"Make a random complex number uniformly distributed in [-1,1]x[-1,1], in the current default precision");
	def("complex_unit", bertini::multiprecision::rand_unit,"Make a random complex number of magnitude 1, in the current default precision");

	
	mpfr_complex (*RandRealNoArgs)() = &bertini::multiprecision::RandomReal;
	def("real_as_complex", RandRealNoArgs, "Make a random real number in [-1,1], as a complex number with imaginary part 0, in the current default precision");
}







} //namespace python
} // namespace bertini