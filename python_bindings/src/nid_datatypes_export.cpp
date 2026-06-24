// python/nid_datatypes_export.cpp — NID result data type registrations (Slice, WitnessSet, NIDResult).
// These are pure data types; no algorithm class_<> instantiations needed.

#include "numerical_irreducible_decomposition_export.hpp"
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/stl_iterator.hpp>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace bertini{
	namespace python{

		using dbl = std::complex<double>;

		// pickle via boost serialization -- so Slice / WitnessSet / NID result round-trip through
		// pickle, copy, and deepcopy (and, in C++, through MPI / threads).  Mirrors the suites in
		// system_export.cpp and mpfr_export.cpp.  WitnessSet restores its system's differentiated
		// state in its own serialization load(), so setstate needs no extra fixup here.
		template<typename T>
		struct BoostSerializePickle : boost::python::pickle_suite
		{
			static boost::python::tuple getinitargs(T const&){ return boost::python::make_tuple(); }

			static boost::python::object getstate(T const& obj)
			{
				std::ostringstream oss;
				{ boost::archive::text_oarchive oa(oss); oa << obj; }
				return boost::python::str(oss.str());
			}

			static void setstate(T& obj, boost::python::object state)
			{
				std::string s = boost::python::extract<std::string>(state)();
				std::istringstream iss(s);
				{ boost::archive::text_iarchive ia(iss); ia >> obj; }
			}
		};

		// ---- Slice: the linear part of a witness set -----------------------------------------
		// A Slice is a stack of linear forms M [x ; 1] (one row per form, the trailing column the
		// constant term), backed by a LinearFormsBlock.  Read-mostly: the factories build one,
		// then it is carried around, evaluated, row-subset, and dropped into a System.
		void ExportSlice(){
			class_<Slice>("Slice", "A linear slice of affine/projective space: a stack of linear forms M [x ; 1].", init<>())
			.def("from_coefficients",
				+[](VariableGroup const& v, bertini::Mat<mpfr_complex> const& aug, bool homogeneous){
					return Slice::FromCoefficients(v, aug, homogeneous);
				},
				(arg("variables"), arg("coefficients"), arg("homogeneous")=false),
				"Build a slice from an augmented coefficient matrix: one row per linear form, num_variables+1 columns, the trailing column carrying each form's constant term (zero, for a homogeneous slice).")
			.staticmethod("from_coefficients")
			.def("random_complex",
				+[](VariableGroup const& v, unsigned dim, bool homogeneous, bool orthogonal){
					return Slice::RandomComplex(v, dim, homogeneous, orthogonal);
				},
				(arg("variables"), arg("dim"), arg("homogeneous")=false, arg("orthogonal")=true),
				"A random complex linear slice of `dim` dimensions (forms) on the given variables.")
			.staticmethod("random_complex")
			.def("random_real",
				+[](VariableGroup const& v, unsigned dim, bool homogeneous, bool orthogonal){
					return Slice::RandomReal(v, dim, homogeneous, orthogonal);
				},
				(arg("variables"), arg("dim"), arg("homogeneous")=false, arg("orthogonal")=true),
				"A random real linear slice of `dim` dimensions (forms) on the given variables.")
			.staticmethod("random_real")
			.def("coefficients",
				+[](Slice const& s){ return bertini::Mat<mpfr_complex>(s.Coefficients()); },
				(arg("self")),
				"The augmented coefficient matrix (one row per linear form, num_variables+1 columns; the trailing column is the constant term).  These rows are ready-made factors for a products-of-linears block.")
			.def("dimension", &Slice::Dimension, (arg("self")), "the dimension of the slice -- the number of linear forms")
			.def("num_variables", &Slice::NumVariables, (arg("self")), "the number of variables the slice is a function of")
			.def("is_homogeneous", &Slice::IsHomogeneous, (arg("self")), "whether the slice was authored without constant terms")
			.def("eval",
				+[](Slice const& s, bertini::Vec<dbl> const& x){ return s.Eval(x); },
				(arg("self"), arg("x")), "evaluate the linear forms at x, in double precision")
			.def("eval",
				+[](Slice const& s, bertini::Vec<mpfr_complex> const& x){ return s.Eval(x); },
				(arg("self"), arg("x")), "evaluate the linear forms at x, in multiple precision")
			.def("add_to", &Slice::AddTo, (arg("self"), arg("system")), "add this slice's linear forms to a System as a linear-forms block")
			.def("head", &Slice::Head, (arg("self"), arg("m")), "a new slice over the same variables, built from the first m linear forms")
			.def("tail", &Slice::Tail, (arg("self"), arg("m")), "a new slice over the same variables, built from the last m linear forms")
			.def("rows",
				+[](Slice const& s, boost::python::list const& indices){
					std::vector<unsigned> idx{
						boost::python::stl_input_iterator<unsigned>(indices),
						boost::python::stl_input_iterator<unsigned>() };
					return s.Rows(idx);
				},
				(arg("self"), arg("indices")), "a new slice over the same variables, built from the chosen linear forms")
			.def("precision", +[](Slice const& s){ return s.Precision(); }, (arg("self")), "get the current working precision of the slice")
			.def("precision", +[](Slice const& s, unsigned p){ s.Precision(p); }, (arg("self"), arg("precision")), "set the working precision of the slice")
			.def_pickle(BoostSerializePickle<Slice>())
			;
		}

		// ---- WitnessSet ----------------------------------------------------------------------
		template<typename NumT>
		void ExportWitnessSet(std::string const& class_name){
			using WS = nag_datatype::WitnessSet<NumT>;
			using PointContT = nag_datatype::PointCont<bertini::Vec<NumT>>;

			class_<WS>(class_name.c_str(), init<>())
			// build all-at-once from points + slice + system, or incrementally (default ctor + add_point/set_*).
			.def("__init__", make_constructor(
				+[](boost::python::list const& pts, Slice const& slc, bertini::System const& sys){
					PointContT points{
						boost::python::stl_input_iterator<bertini::Vec<NumT>>(pts),
						boost::python::stl_input_iterator<bertini::Vec<NumT>>() };
					return new WS(points, slc, sys);
				}),
				"construct a witness set from a list of witness points, a slice, and a system")
			.def("degree", &WS::Degree, "the degree of the component, i.e. the number of witness points")
			.def("dimension", &WS::Dimension, "the dimension of the component")
			.def("is_consistent", &WS::IsConsistent, "whether the slice dimension matches the system's underdeterminedness")
			.def("get_point", &WS::GetPoint, return_value_policy<copy_const_reference>(), "get the i-th witness point")
			.def("get_points",
				+[](WS const& w){
					boost::python::list out;
					for (auto const& p : w.GetPoints())
						out.append(p);
					return out;
				},
				"the witness points, as a list")
			.def("get_slice", &WS::GetSlice, return_internal_reference<>(), "the linear slice of this witness set")
			.def("get_system", &WS::GetSystem, return_internal_reference<>(), "the system this witness set is for")
			// mutators -- build a witness set as you go
			.def("add_point", &WS::AddPoint, (arg("self"), arg("point")), "append a witness point")
			.def("set_points",
				+[](WS& w, boost::python::list const& pts){
					PointContT points{
						boost::python::stl_input_iterator<bertini::Vec<NumT>>(pts),
						boost::python::stl_input_iterator<bertini::Vec<NumT>>() };
					w.SetPoints(points);
				},
				(arg("self"), arg("points")), "replace the witness points with the given list")
			.def("set_slice", &WS::SetSlice, (arg("self"), arg("slice")), "set the linear slice")
			.def("set_system", &WS::SetSystem, (arg("self"), arg("system")), "set the system")
			.def_pickle(BoostSerializePickle<WS>())
			;
		}

		template<typename NumT>
		void ExportNIDResult(std::string const& class_name){
			using R = nag_datatype::NumericalIrreducibleDecomposition<NumT>;
			class_<R>(class_name.c_str(), init<>())
			.def("nonempty_codimensions", &R::NonEmptyCodimensions, "the distinct codimensions which contain at least one component")
			.def("num_witness_sets", &R::NumWitnessSets, "the number of stored witness sets")
			.def("get_witness_set", &R::GetWitnessSet, return_internal_reference<>(), "get the i-th stored witness set")
			.def_pickle(BoostSerializePickle<R>())
			;
		}

		void ExportNIDDataTypes(){
			ExportSlice();
			ExportWitnessSet<dbl_complex>("WitnessSetDoublePrecision");
			ExportWitnessSet<mpfr_complex>("WitnessSetMultiplePrecision");
			ExportNIDResult<dbl_complex>("NumericalIrreducibleDecompositionDoublePrecision");
			ExportNIDResult<mpfr_complex>("NumericalIrreducibleDecompositionMultiplePrecision");
		}

}} // namespaces
