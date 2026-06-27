#include "random_export.hpp"
#include "bertini2/random.hpp"


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

	complex_mp (*RandRealNoArgs)() = &bertini::multiprecision::RandomReal;
	def("real_as_complex", RandRealNoArgs, "Make a random real number in [-1,1], as a complex number with imaginary part 0, in the current default precision");

	def("set_random_seed", &bertini::SetGlobalSeed, boost::python::arg("seed") = 0ul,
		"Set the global RNG seed (0 = draw from entropy). Call before constructing any "
		"homotopy or solver to get reproducible results. The effective seed (which may "
		"differ from 0 when entropy is used) is retrievable via get_random_seed().");
	def("get_random_seed", &bertini::GetGlobalSeed,
		"Return the effective global RNG seed. If set_random_seed has not been called, "
		"draws from entropy on first call and caches the result.");
}







} //namespace python
} // namespace bertini