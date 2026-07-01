#include "random_export.hpp"
#include "bertini2/random.hpp"
#include "bertini2/num_traits.hpp"
#include "bertini2/eigen_extensions.hpp"


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

	complex_mp (*RandComplexBM)() = &bertini::multiprecision::RandomComplexBoundedModulus;
	def("complex_bounded_modulus", RandComplexBM, "Make a random complex number of bounded modulus (box-uniform in [-1,1]x[-1,1], magnitude at most sqrt(2), away from 0 -- the Bertini coefficient draw), in the current default precision");

	complex_mp (*RandRealBM)() = &bertini::multiprecision::RandomRealBoundedModulus;
	def("real_bounded_modulus", RandRealBM, "Make a random real number of bounded modulus (box-uniform in [-1,1]), as a complex number with imaginary part 0, in the current default precision");

	def("real_unit", +[]() { return complex_mp(bertini::RandomUnit<bertini::real_mp>()); },
		"Make a random real number of unit modulus (i.e. +1 or -1), as a complex number with imaginary part 0, in the current default precision");

	def("conjugate_orthonormal_matrix",
		+[](unsigned rows, unsigned cols, bool real) -> bertini::Mat<complex_mp> {
			if (real) {
				bertini::Mat<bertini::real_mp> M = bertini::RandomConjugateOrthonormalMatrix<bertini::real_mp>(rows, cols);
				return M.unaryExpr([](bertini::real_mp const& r){ return complex_mp(r); });
			}
			return bertini::RandomConjugateOrthonormalMatrix<complex_mp>(rows, cols);
		},
		(arg("rows"), arg("cols"), arg("real")=false),
		"A rows x cols random conjugate-orthonormal matrix (its rows orthonormal -- QR-factored from a square matrix of units, then truncated; perfectly conditioned), at the current default precision.  real=True yields a real orthogonal matrix (entries with zero imaginary part).  Returns a complex_mp matrix.");

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