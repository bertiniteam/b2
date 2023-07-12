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



	def("complex_in_minus_one_to_one", bertini::multiprecision::rand);

}







} //namespace python
} // namespace bertini