


#pragma once


#ifndef BERTINI_PYTHON_EIGENPY_INTERACTION_HPP
#define BERTINI_PYTHON_EIGENPY_INTERACTION_HPP


#include "python_common.hpp"

#include <eigenpy/eigenpy.hpp>
#include <eigenpy/user-type.hpp>
#include <eigenpy/ufunc.hpp>





// this code derived from 
// https://github.com/stack-of-tasks/eigenpy/issues/365
// where I asked about using custom types, and @jcarpent responded with a discussion
// of an application of this in Pinnochio, a library for rigid body dynamics.
namespace eigenpy
{
  namespace internal
  {


// template specialization for real numbers
template <>
struct getitem<bertini::mpfr_float>
{
 using NumT = bertini::mpfr_float;

  static PyObject* run(void* data, void* /* arr */) {
    NumT & mpfr_scalar = *static_cast<NumT*>(data);
    auto & backend = mpfr_scalar.backend();
   
    if(backend.data()[0]._mpfr_d == 0) // If the mpfr_scalar is not initialized, we have to init it.  
    {
      mpfr_scalar = NumT(0);
    }
    boost::python::object m(boost::ref(mpfr_scalar));
    Py_INCREF(m.ptr());
    return m.ptr();
  }
};




// a template specialization for complex numbers
  template <>
  struct getitem<bertini::mpfr_complex>
  {
   using NumT = bertini::mpfr_complex;

    static PyObject* run(void* data, void* /* arr */) {
      NumT & mpfr_scalar = *static_cast<NumT*>(data);
      auto & backend = mpfr_scalar.backend();
     
      if(backend.data()[0].re->_mpfr_d == 0) // If the mpfr_scalar is not initialized, we have to init it.  
      {
        mpfr_scalar = NumT(0);
      }
      boost::python::object m(boost::ref(mpfr_scalar));
      Py_INCREF(m.ptr());
      return m.ptr();
    }
  };


} // namespace internal



// i lifted this from EigenPy and adapted it, basically removing the calls for the comparitors.
template <typename Scalar>
void registerUfunct_without_comparitors(){
	  const int type_code = Register::getTypeCode<Scalar>();

	  PyObject *numpy_str;
	#if PY_MAJOR_VERSION >= 3
	  numpy_str = PyUnicode_FromString("numpy");
	#else
	  numpy_str = PyString_FromString("numpy");
	#endif
	  PyObject *numpy;
	  numpy = PyImport_Import(numpy_str);
	  Py_DECREF(numpy_str);

	  import_ufunc();

	  // Matrix multiply
	  {
	    int types[3] = {type_code, type_code, type_code};

	    std::stringstream ss;
	    ss << "return result of multiplying two matrices of ";
	    ss << bp::type_info(typeid(Scalar)).name();
	    PyUFuncObject *ufunc =
	        (PyUFuncObject *)PyObject_GetAttrString(numpy, "matmul");
	    if (!ufunc) {
	      std::stringstream ss;
	      ss << "Impossible to define matrix_multiply for given type "
	         << bp::type_info(typeid(Scalar)).name() << std::endl;
	      eigenpy::Exception(ss.str());
	    }
	    if (PyUFunc_RegisterLoopForType((PyUFuncObject *)ufunc, type_code,
	                                    &internal::gufunc_matrix_multiply<Scalar>,
	                                    types, 0) < 0) {
	      std::stringstream ss;
	      ss << "Impossible to register matrix_multiply for given type "
	         << bp::type_info(typeid(Scalar)).name() << std::endl;
	      eigenpy::Exception(ss.str());
	    }

	    Py_DECREF(ufunc);
	  }

	  // Binary operators
	  EIGENPY_REGISTER_BINARY_UFUNC(add, type_code, Scalar, Scalar, Scalar);
	  EIGENPY_REGISTER_BINARY_UFUNC(subtract, type_code, Scalar, Scalar, Scalar);
	  EIGENPY_REGISTER_BINARY_UFUNC(multiply, type_code, Scalar, Scalar, Scalar);
	  EIGENPY_REGISTER_BINARY_UFUNC(divide, type_code, Scalar, Scalar, Scalar);

	  // Comparison operators
	  EIGENPY_REGISTER_BINARY_UFUNC(equal, type_code, Scalar, Scalar, bool);
	  EIGENPY_REGISTER_BINARY_UFUNC(not_equal, type_code, Scalar, Scalar, bool);

	  //these are commented out because the comparisons are NOT defined for complex types!!
	  // EIGENPY_REGISTER_BINARY_UFUNC(greater, type_code, Scalar, Scalar, bool);
	  // EIGENPY_REGISTER_BINARY_UFUNC(less, type_code, Scalar, Scalar, bool);
	  // EIGENPY_REGISTER_BINARY_UFUNC(greater_equal, type_code, Scalar, Scalar, bool);
	  // EIGENPY_REGISTER_BINARY_UFUNC(less_equal, type_code, Scalar, Scalar, bool);

	  // Unary operators
	  EIGENPY_REGISTER_UNARY_UFUNC(negative, type_code, Scalar, Scalar);

	  Py_DECREF(numpy);
}

} // namespace eigenpy



namespace bertini{
	namespace python{

void ExportEigenPy();
	
}} // namespaces


#endif // include guard

