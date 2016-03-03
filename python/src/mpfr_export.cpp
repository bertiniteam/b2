




#include "mpfr_export.hpp"



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
			
			.def("__str__", &MPFRBaseVisitor::__str__)
			;
			
			
			def("abs", &MPFRBaseVisitor::__abs__);
			
		}
		
		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRFloatBaseVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def("__div__",&MPFRFloatBaseVisitor::__div__).def("__idiv__",&MPFRFloatBaseVisitor::__idiv__)
			.def("__pow__",&MPFRFloatBaseVisitor::__pow__)
			
			.def("precision", bmpprec1)
			.def("precision", bmpprec2)
			;
			
			def("exp", &MPFRFloatBaseVisitor::__exp__);
			def("log", &MPFRFloatBaseVisitor::__log__);
			def("sqrt", &MPFRFloatBaseVisitor::__sqrt__);
			
			def("sin", &MPFRFloatBaseVisitor::__sin__);
			def("cos", &MPFRFloatBaseVisitor::__cos__);
			def("tan", &MPFRFloatBaseVisitor::__tan__);
			
			def("asin", &MPFRFloatBaseVisitor::__asin__);
			def("acos", &MPFRFloatBaseVisitor::__acos__);
			def("atan", &MPFRFloatBaseVisitor::__atan__);
			
			def("sinh", &MPFRFloatBaseVisitor::__sinh__);
			def("cosh", &MPFRFloatBaseVisitor::__cosh__);
			def("tanh", &MPFRFloatBaseVisitor::__tanh__);
			
		}

		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRIntVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def("__add__",&MPFRIntVisitor::__add_int).def("__iadd__",&MPFRIntVisitor::__iadd_int)
			.def("__radd__",&MPFRIntVisitor::__radd_int)
			.def("__sub__",&MPFRIntVisitor::__sub_int).def("__isub__",&MPFRIntVisitor::__isub_int)
			.def("__rsub__",&MPFRIntVisitor::__rsub_int)
			.def("__mul__",&MPFRIntVisitor::__mul_int).def("__imul__",&MPFRIntVisitor::__imul_int)
			.def("__rmul__",&MPFRIntVisitor::__rmul_int)
			.def("__pow__",&MPFRIntVisitor::__pow__)
			.def("__repr__", &MPFRIntVisitor::__repr__)

			.def(self < self)
			.def(self <= self)
			.def(self > self)
			.def(self >= self)
			.def(self == self)
			.def(self != self)
			;

		};

		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRRationalVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def("__add__",&MPFRRationalVisitor::__add_int).def("__iadd__",&MPFRRationalVisitor::__iadd_int)
			.def("__radd__",&MPFRRationalVisitor::__radd_int)
			.def("__sub__",&MPFRRationalVisitor::__sub_int).def("__isub__",&MPFRRationalVisitor::__isub_int)
			.def("__rsub__",&MPFRRationalVisitor::__rsub_int)
			.def("__mul__",&MPFRRationalVisitor::__mul_int).def("__imul__",&MPFRRationalVisitor::__imul_int)
			.def("__rmul__",&MPFRRationalVisitor::__rmul_int)
			.def("__div__",&MPFRRationalVisitor::__div__).def("__idiv__",&MPFRRationalVisitor::__idiv__)
			.def("__repr__", &MPFRRationalVisitor::__repr__)
			
			.def(self < self)
			.def(self <= self)
			.def(self > self)
			.def(self >= self)
			.def(self == self)
			.def(self != self)
			;
			
		};

		
		
		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRFloatVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRFloatBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
			.def("__add__",&MPFRFloatVisitor::__add_int).def("__iadd__",&MPFRFloatVisitor::__iadd_int)
			.def("__radd__",&MPFRFloatVisitor::__radd_int)
			.def("__sub__",&MPFRFloatVisitor::__sub_int).def("__isub__",&MPFRFloatVisitor::__isub_int)
			.def("__rsub__",&MPFRFloatVisitor::__rsub_int)
			.def("__mul__",&MPFRFloatVisitor::__mul_int).def("__imul__",&MPFRFloatVisitor::__imul_int)
			.def("__rmul__",&MPFRFloatVisitor::__rmul_int)
			.def("__div__",&MPFRFloatVisitor::__div_int).def("__idiv__",&MPFRFloatVisitor::__idiv_int)
			.def("__rdiv__",&MPFRFloatVisitor::__rdiv_int)
			.def("__pow__",&MPFRFloatVisitor::__pow_int)
			
			.def("__repr__", &MPFRFloatVisitor::__repr__)
			
			.def(self < self)
			.def(self <= self)
			.def(self > self)
			.def(self >= self)
			.def(self == self)
			.def(self != self)
			;
			
			def("default_precision", MPFRFloatVisitor<bmp>::def_prec1);
			def("default_precision", MPFRFloatVisitor<bmp>::def_prec2);

		};


		template<typename MPFRBaseT>
		template<class PyClass>
		void MPFRComplexVisitor<MPFRBaseT>::visit(PyClass& cl) const
		{
			MPFRFloatBaseVisitor<MPFRBaseT>().visit(cl);
			
			cl
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
			class_<boost::multiprecision::mpz_int>("mpfr_int", init<>())
			.def(init<int>())
			.def(init<boost::multiprecision::mpz_int>())
			.def(MPFRIntVisitor<boost::multiprecision::mpz_int>())
			;

			class_<boost::multiprecision::mpq_rational>("mpfr_rational", init<>())
			.def(init<int, int>())
			.def(init<boost::multiprecision::mpz_int,boost::multiprecision::mpz_int>())
			.def(init<boost::multiprecision::mpq_rational>())
			.def(MPFRRationalVisitor<boost::multiprecision::mpq_rational>())
			;

			class_<bmp>("mpfr_float", init<>())
			.def(init<std::string>())
			.def(init<int>())
			.def(init<long int>())
			.def(init<bmp>())
			.def(MPFRFloatVisitor<bmp>())
			;
			
			class_<bertini::complex>("mpfr_complex", init<>())
			.def(init<double>())
			.def(init<mpfr_float>())
			.def(init<std::string>())
			.def(init<mpfr_float,mpfr_float>())
			.def(init<double, double>())
			.def(init<std::string, mpfr_float>())
			.def(init<mpfr_float, std::string>())
			.def(init<std::string, std::string>())
			.def(init<bertini::complex>())
			.def(MPFRComplexVisitor<bertini::complex>())
			;
			
			
			
		};

		
	} //namespace python
} // namespace bertini


