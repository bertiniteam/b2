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

	// --- friendly factories (issue #294): make a random real / complex / vector directly ---
	// Continuous bounded-modulus draws (away from 0 and infinity -- the Bertini genericity draw),
	// so a random projection / linear functional is generic and changes with set_random_seed.  Unlike
	// the orthonormal random_matrix, these are NOT quantized: a real random_vector is the right tool
	// for a generic real projection direction.
	def("random_real",
		+[]() -> bertini::real_mp { return bertini::multiprecision::RandomRealBoundedModulus().real(); },
		"Make a random real number (real_mp) of bounded modulus (box-uniform in [-1,1], away from 0), at the current default precision.  Reproducible via set_random_seed.");

	def("random_complex",
		+[]() -> complex_mp { return bertini::multiprecision::RandomComplexBoundedModulus(); },
		"Make a random complex number (complex_mp) of bounded modulus (away from 0 and infinity), at the current default precision.  Reproducible via set_random_seed.");

	def("random_vector",
		+[](unsigned size, bool real) -> object {
			if (real) {
				bertini::Vec<bertini::real_mp> v(size);
				for (unsigned i = 0; i < size; ++i)
					v(i) = bertini::multiprecision::RandomRealBoundedModulus().real();
				return object(v);
			}
			bertini::Vec<complex_mp> v(size);
			for (unsigned i = 0; i < size; ++i)
				v(i) = bertini::multiprecision::RandomComplexBoundedModulus();
			return object(v);
		},
		(arg("size"), arg("real") = false),
		"A random length-`size` vector of bounded-modulus numbers -- real_mp when real=True, else complex_mp -- at the current default precision.  The natural random projection / linear-functional coefficient vector (generic and seed-reproducible, unlike the quantized orthonormal random_matrix).");

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
	def("derive_solve_seed", &bertini::DeriveSolveSeed,
		"Capture the session stream's current position as an effective per-solve seed: "
		"one draw off the current stream (advancing it), nonzero, 32-bit portable. "
		"bertini.solve uses this when no explicit seed is given, so a recorded run's "
		"seed reproduces that run standalone -- deterministic from the session master "
		"when one was set, random for a never-seeded session.");
}







} //namespace python
} // namespace bertini