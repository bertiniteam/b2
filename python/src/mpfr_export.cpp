




#include "mpfr_visitors.hpp"



namespace bertini{
	namespace python{
		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRBaseVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			cl
			.def("__neg__",&MPFRBaseVisitor::__neg__)
			
			.def("__add__",&MPFRBaseVisitor::__add__).def("__iadd__",&MPFRBaseVisitor::__iadd__)
			.def("__sub__",&MPFRBaseVisitor::__sub__).def("__isub__",&MPFRBaseVisitor::__isub__)
			.def("__mul__",&MPFRBaseVisitor::__mul__).def("__imul__",&MPFRBaseVisitor::__imul__)
			.def("__div__",&MPFRBaseVisitor::__div__).def("__idiv__",&MPFRBaseVisitor::__idiv__)
			.def("__pow__",&MPFRBaseVisitor::__pow__)
			
			.def("precision", bmpprec1)
			.def("precision", bmpprec2)
			;
			
			
			def("abs", &MPFRBaseVisitor::__abs__);
			def("exp", &MPFRBaseVisitor::__exp__);
			def("log", &MPFRBaseVisitor::__log__);
			def("sqrt", &MPFRBaseVisitor::__sqrt__);
			
			def("sin", &MPFRBaseVisitor::__sin__);
			def("cos", &MPFRBaseVisitor::__cos__);
			def("tan", &MPFRBaseVisitor::__tan__);
			
			def("asin", &MPFRBaseVisitor::__asin__);
			def("acos", &MPFRBaseVisitor::__acos__);
			def("atan", &MPFRBaseVisitor::__atan__);
			
			def("sinh", &MPFRBaseVisitor::__sinh__);
			def("cosh", &MPFRBaseVisitor::__cosh__);
			def("tanh", &MPFRBaseVisitor::__tanh__);
			
		}

		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRFloatVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def(init<std::string>())
			.def(init<int>())
			.def(init<long int>())
			.def(init<MPFRBaseT>())
			
			.def("__add__",&MPFRFloatVisitor::__add_int).def("__iadd__",&MPFRFloatVisitor::__iadd_int)
			.def("__radd__",&MPFRFloatVisitor::__radd_int)
			.def("__sub__",&MPFRFloatVisitor::__sub_int).def("__isub__",&MPFRFloatVisitor::__isub_int)
			.def("__rsub__",&MPFRFloatVisitor::__rsub_int)
			.def("__mul__",&MPFRFloatVisitor::__mul_int).def("__imul__",&MPFRFloatVisitor::__imul_int)
			.def("__rmul__",&MPFRFloatVisitor::__rmul_int)
			.def("__div__",&MPFRFloatVisitor::__div_int).def("__idiv__",&MPFRFloatVisitor::__idiv_int)
			.def("__rdiv__",&MPFRFloatVisitor::__rdiv_int)
			.def("__pow__",&MPFRFloatVisitor::__pow_int)
			
			.def("__str__", &MPFRFloatVisitor::__str__)
			.def("__repr__", &MPFRFloatVisitor::__repr__)
			
			.def(self < self)
			.def(self <= self)
			.def(self > self)
			.def(self >= self)
			.def(self == self)
			.def(self != self)
			;
			
			def("default_precision", def_prec1);
			def("default_precision", def_prec2);
		};


		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRComplexVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def(init<double>())
			.def(init<mpfr_float>())
			.def(init<std::string>())
			.def(init<mpfr_float,mpfr_float>())
			.def(init<double, double>())
			.def(init<std::string, mpfr_float>())
			.def(init<mpfr_float, std::string>())
			.def(init<std::string, std::string>())
			.def(init<MPFRBaseT>())
			
			.def("__add__",&MPFRComplexVisitor::__add_float).def("__iadd__",&MPFRComplexVisitor::__iadd_float)
			.def("__radd__",&MPFRComplexVisitor::__radd_float)
			.def("__sub__",&MPFRComplexVisitor::__sub_float).def("__isub__",&MPFRComplexVisitor::__isub_float)
			.def("__rsub__",&MPFRComplexVisitor::__rsub_float)
			.def("__mul__",&MPFRComplexVisitor::__mul_float).def("__imul__",&MPFRComplexVisitor::__imul_float)
			.def("__rmul__",&MPFRComplexVisitor::__rmul_float)
			.def("__div__",&MPFRComplexVisitor::__div_float).def("__idiv__",&MPFRComplexVisitor::__idiv_float)
			.def("__rdiv__",&MPFRComplexVisitor::__rdiv_float)
			.def("__pow__",&MPFRComplexVisitor::__pow_int)
			.def("__pow__",&MPFRComplexVisitor::__pow_float)
			
			.def("__str__", &MPFRComplexVisitor::__str__)
			.def("__repr__", &MPFRComplexVisitor::__repr__)
			
			.add_property("real", getreal, setreal)
			.add_property("imag", getimag, setimag)
			;
			
			
			def("real",&real);
			def("imag",&imag);
			
			def("abs2",&MPFRBaseT::abs2);
			def("polar",&polar);
			def("norm",&MPFRBaseT::norm);
			def("conj",&MPFRBaseT::conj);
			def("arg",&arg);
			
			def("square",&square);
			def("cube",&cube);
			def("inverse", &inverse);
			def("asinh",&asinh);
			def("acosh",&acosh);
			def("atanh",&atanh);
			
		}


		void ExportMpfr()
		{
			class_<bmp>("mpfr_float", init<>())
			.def(MPFRFloatVisitor<bmp>());
			
			class_<bertini::complex>("mpfr_complex", init<>())
			.def(MPFRComplexVisitor<bertini::complex>());
		};

		
	} //namespace python
} // namespace bertini


