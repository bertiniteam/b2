

#pragma once

#ifndef BERTINI_PYTHON_EIGENPY_INTERACTION_HPP
#define BERTINI_PYTHON_EIGENPY_INTERACTION_HPP

#include "python_common.hpp"

#include <cstring>

#include <eigenpy/eigenpy.hpp>
#include <eigenpy/user-type.hpp>
#include <eigenpy/ufunc.hpp>

// this code derived from
// https://github.com/stack-of-tasks/eigenpy/issues/365
// where I asked about using custom types, and @jcarpent responded with a discussion
// of an application of this in Pinnochio, a library for rigid body dynamics.
//
// ---------------------------------------------------------------------------
// On uninitialized array slots (the cause of intermittent SIGABRT/SIGSEGV):
//
// numpy zero-fills freshly allocated buffers for these dtypes (eigenpy sets
// NPY_NEEDS_INIT when registering, since eigenpy 2.6.4).  An all-zero
// mpfr_t/mpc_t is NOT a valid value — it is Boost.Multiprecision's
// "uninitialized" sentinel (_mpfr_d == 0).  BMP's own assignment operators
// check the sentinel and initialize first, so WRITES into fresh slots are
// safe.  But anything that READS a never-written slot and hands it straight
// to mpfr/mpc crashes inside libmpfr/libmpc.  Three defenses live here:
//
//   1. getitem specializations heal a zeroed slot in place on read.
//   2. zeroinit_setitem/HardenSetitem zero the slot before delegating, so
//      malloc-dirty memory (seen in manylinux, ADR-0003) cannot defeat BMP's
//      null check on the write path.
//   3. guarded ufunc loops and cast specializations substitute an exact zero
//      on the read side (value_or_zero), replacing eigenpy's unguarded loops
//      which segfault on np.zeros/np.empty slots that were never written.
// ---------------------------------------------------------------------------
namespace eigenpy
{
	namespace internal
	{

		// trait: detect Boost.Multiprecision's "uninitialized" sentinel in a
		// zero-filled numpy slot.
		template <typename T>
		struct mpfr_slot
		{
			static bool uninitialized(T const&) { return false; }
		};

		template <>
		struct mpfr_slot<bertini::real_mp>
		{
			static bool uninitialized(bertini::real_mp const& x)
			{
				return x.backend().data()[0]._mpfr_d == 0;
			}
		};

		template <>
		struct mpfr_slot<bertini::complex_mp>
		{
			static bool uninitialized(bertini::complex_mp const& x)
			{
				return x.backend().data()[0].re->_mpfr_d == 0;
			}
		};

		// read-side guard: yield `zero` for an uninitialized slot, the slot's value otherwise.
		template <typename T>
		inline T const& value_or_zero(T const& x, T const& zero)
		{
			return mpfr_slot<T>::uninitialized(x) ? zero : x;
		}


		// template specialization for real numbers
		template <>
		struct getitem<bertini::real_mp>
		{
			using NumT = bertini::real_mp;

			static PyObject *run(void *data, void * /* arr */)
			{
				NumT &mpfr_scalar = *static_cast<NumT *>(data);

				if (mpfr_slot<NumT>::uninitialized(mpfr_scalar)) // heal a never-written slot in place
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
		struct getitem<bertini::complex_mp>
		{
			using NumT = bertini::complex_mp;

			static PyObject *run(void *data, void * /* arr */)
			{
				NumT &mpfr_scalar = *static_cast<NumT *>(data);

				if (mpfr_slot<NumT>::uninitialized(mpfr_scalar)) // heal a never-written slot in place
				{
					mpfr_scalar = NumT(0);
				}
				boost::python::object m(boost::ref(mpfr_scalar));
				Py_INCREF(m.ptr());
				return m.ptr();
			}
		};

		// Zero-initialization guard for numpy's setitem on MPFR-backed dtypes.
		//
		// eigenpy's SpecialMethods<T>::setitem copy-assigns into the raw numpy slot:
		//     T& dest = *static_cast<T*>(dest_ptr);
		//     dest = src;
		// Boost.Multiprecision's operator= decides whether the destination needs
		// mpfr/mpc initialization by testing _mpfr_d == nullptr (see
		// mpc_complex_imp::operator= in boost/multiprecision/mpc.hpp).  numpy is
		// supposed to hand us zeroed memory — eigenpy registers the dtype with
		// NPY_NEEDS_INIT — but the manylinux containers have been observed to
		// deliver malloc-dirty slots (ADR-0003): garbage non-null _mpfr_d defeats
		// the null check, mpc_set runs on garbage, and MPFR_ASSERTN aborts.
		//
		// Zeroing the slot before delegating forces operator= onto its
		// init-before-set path regardless of what the allocator delivered.  If the
		// slot held a previously initialized value, its mpfr allocation leaks; numpy
		// never destructs user-dtype elements anyway (every discarded array of this
		// dtype already leaks its elements), so this converts a crash on dirty
		// memory into a small, bounded leak on element overwrite.
		//
		// This is a wrapper installed post-registration (see HardenSetitem below)
		// rather than a template specialization because, unlike getitem,
		// eigenpy provides no customization point for setitem.
		template <typename NumT>
		struct zeroinit_setitem
		{
			static inline PyArray_SetItemFunc *original = nullptr;

			static int run(PyObject *src_obj, void *dest_ptr, void *array)
			{
				std::memset(dest_ptr, 0, sizeof(NumT));
				return original(src_obj, dest_ptr, array);
			}
		};


		// ----- guarded ufunc loops ------------------------------------------------
		// These replace eigenpy's EIGENPY_REGISTER_{BINARY,UNARY}_UFUNC loop bodies
		// (and its gufunc_matrix_multiply), which read input slots unguarded and
		// segfault inside libmpfr/libmpc on never-written np.zeros/np.empty slots.
		// Writes into the output slot go through BMP operator=, which initializes
		// a zeroed destination itself.

		struct op_add           { template <typename T> static T    apply(T const& x, T const& y) { return T(x + y); } };
		struct op_subtract      { template <typename T> static T    apply(T const& x, T const& y) { return T(x - y); } };
		struct op_multiply      { template <typename T> static T    apply(T const& x, T const& y) { return T(x * y); } };
		struct op_divide        { template <typename T> static T    apply(T const& x, T const& y) { return T(x / y); } };
		struct op_equal         { template <typename T> static bool apply(T const& x, T const& y) { return x == y; } };
		struct op_not_equal     { template <typename T> static bool apply(T const& x, T const& y) { return x != y; } };
		struct op_greater       { template <typename T> static bool apply(T const& x, T const& y) { return x > y; } };
		struct op_less          { template <typename T> static bool apply(T const& x, T const& y) { return x < y; } };
		struct op_greater_equal { template <typename T> static bool apply(T const& x, T const& y) { return x >= y; } };
		struct op_less_equal    { template <typename T> static bool apply(T const& x, T const& y) { return x <= y; } };

		struct op_negative { template <typename T> static T apply(T const& x) { return T(-x); } };
		struct op_square   { template <typename T> static T apply(T const& x) { return T(x * x); } };
		struct op_sqrt
		{
			template <typename T> static T apply(T const& x)
			{
				using boost::multiprecision::sqrt;
				using std::sqrt;
				return T(sqrt(x));
			}
		};

		template <typename T, typename Op>
		void guarded_binary_op(
				char **args, EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *dimensions,
				EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *steps, void * /*data*/)
		{
			npy_intp is0 = steps[0], is1 = steps[1], os = steps[2], n = *dimensions;
			char *i0 = args[0], *i1 = args[1], *o = args[2];
			const T zero(0);
			for (npy_intp k = 0; k < n; ++k)
			{
				T const& x = value_or_zero(*reinterpret_cast<T const*>(i0), zero);
				T const& y = value_or_zero(*reinterpret_cast<T const*>(i1), zero);
				T& res = *reinterpret_cast<T*>(o);
				res = Op::apply(x, y);
				i0 += is0;
				i1 += is1;
				o += os;
			}
		}

		template <typename T, typename Op>
		void guarded_compare_op(
				char **args, EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *dimensions,
				EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *steps, void * /*data*/)
		{
			npy_intp is0 = steps[0], is1 = steps[1], os = steps[2], n = *dimensions;
			char *i0 = args[0], *i1 = args[1], *o = args[2];
			const T zero(0);
			for (npy_intp k = 0; k < n; ++k)
			{
				T const& x = value_or_zero(*reinterpret_cast<T const*>(i0), zero);
				T const& y = value_or_zero(*reinterpret_cast<T const*>(i1), zero);
				bool& res = *reinterpret_cast<bool*>(o);
				res = Op::apply(x, y);
				i0 += is0;
				i1 += is1;
				o += os;
			}
		}

		template <typename T, typename Op>
		void guarded_unary_op(
				char **args, EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *dimensions,
				EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *steps, void * /*data*/)
		{
			npy_intp is = steps[0], os = steps[1], n = *dimensions;
			char *i = args[0], *o = args[1];
			const T zero(0);
			for (npy_intp k = 0; k < n; ++k)
			{
				T const& x = value_or_zero(*reinterpret_cast<T const*>(i), zero);
				T& res = *reinterpret_cast<T*>(o);
				res = Op::apply(x);
				i += is;
				o += os;
			}
		}

		// guarded matmul: mirrors eigenpy::internal::{matrix_multiply,gufunc_matrix_multiply}
		// stride logic, with the inner dot product reading through value_or_zero.
		template <typename T>
		void guarded_matrix_multiply(char **args, npy_intp const *dimensions,
		                             npy_intp const *steps)
		{
			char *ip1 = args[0], *ip2 = args[1], *op = args[2];
			npy_intp dm = dimensions[0], dn = dimensions[1], dp = dimensions[2];
			npy_intp is1_m = steps[0], is1_n = steps[1], is2_n = steps[2],
			         is2_p = steps[3], os_m = steps[4], os_p = steps[5];

			const T zero(0);
			for (npy_intp m = 0; m < dm; ++m)
			{
				for (npy_intp p = 0; p < dp; ++p)
				{
					T sum(0);
					char *a = ip1, *b = ip2;
					for (npy_intp k = 0; k < dn; ++k)
					{
						T const& x = value_or_zero(*reinterpret_cast<T const*>(a), zero);
						T const& y = value_or_zero(*reinterpret_cast<T const*>(b), zero);
						sum += x * y;
						a += is1_n;
						b += is2_n;
					}
					T& res = *reinterpret_cast<T*>(op);
					res = sum;
					ip2 += is2_p;
					op += os_p;
				}
				ip2 -= is2_p * dp;
				op -= os_p * dp;
				ip1 += is1_m;
				op += os_m;
			}
		}

		template <typename T>
		void guarded_gufunc_matrix_multiply(
				char **args, EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *dimensions,
				EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *steps, void * /*func*/)
		{
			npy_intp dN = dimensions[0];
			npy_intp s0 = steps[0], s1 = steps[1], s2 = steps[2];
			char *args_local[3] = {args[0], args[1], args[2]};
			for (npy_intp N_ = 0; N_ < dN; ++N_)
			{
				guarded_matrix_multiply<T>(args_local, dimensions + 1, steps + 3);
				args_local[0] += s0;
				args_local[1] += s1;
				args_local[2] += s2;
			}
		}

		// guarded numpy dot/inner (the PyArray_ArrFuncs `dotfunc` slot).
		//
		// eigenpy's SpecialMethods<T>::dotfunc maps both operands as Eigen vectors
		// and calls v0.dot(v1), reading every slot unguarded.  np.dot / np.inner /
		// 1-D '@' on a never-written np.zeros/np.empty array therefore feed an
		// all-zero sentinel (mpfr/mpc _mpfr_d == 0) straight into libmpfr and abort
		// (MPFR_ASSERTN — observed both as init2.c "p>=1" and mpfr_abort_prec_max).
		// numpy zero-fills NEEDS_INIT buffers, so dot operands are always either a
		// valid value or the all-zero sentinel (never malloc-dirty), so reading
		// through value_or_zero is sufficient.  Mirrors guarded_matrix_multiply:
		// accumulate at the ambient default precision and write through operator=
		// into the (zero-filled) scalar output slot.
		template <typename T>
		void guarded_dotfunc(void *ip0_, npy_intp is0, void *ip1_, npy_intp is1,
		                     void *op, npy_intp n, void * /*arr*/)
		{
			const T zero(0);
			T acc(0);
			char *p0 = static_cast<char*>(ip0_);
			char *p1 = static_cast<char*>(ip1_);
			for (npy_intp i = 0; i < n; ++i)
			{
				T const& x = value_or_zero(*reinterpret_cast<T const*>(p0), zero);
				T const& y = value_or_zero(*reinterpret_cast<T const*>(p1), zero);
				acc += x * y;
				p0 += is0;
				p1 += is1;
			}
			*reinterpret_cast<T*>(op) = acc;
		}

	} // namespace internal


	// Guard the registered numpy cast loops too: eigenpy's internal::cast does
	// `to[i] = eigenpy::cast<From,To>::run(from[i])`, reading From slots
	// unguarded.  The primary template is documented as specializable.
	template <typename To>
	struct cast<bertini::real_mp, To>
	{
		static To run(bertini::real_mp const& from)
		{
			if (internal::mpfr_slot<bertini::real_mp>::uninitialized(from))
				return To(0);
			return static_cast<To>(from);
		}
	};

	template <typename To>
	struct cast<bertini::complex_mp, To>
	{
		static To run(bertini::complex_mp const& from)
		{
			if (internal::mpfr_slot<bertini::complex_mp>::uninitialized(from))
				return To(0);
			return static_cast<To>(from);
		}
	};


	// Install the zero-initialization setitem guard for an MPFR-backed dtype.
	// Call immediately after eigenpy::registerNewType<NumT>(), before any arrays
	// of this dtype can exist.  See internal::zeroinit_setitem for the rationale.
	template <typename NumT>
	void HardenSetitem()
	{
		PyArray_Descr *descr = Register::getPyArrayDescr<NumT>();
		PyArray_ArrFuncs *funcs = PyDataType_GetArrFuncs(descr);
		internal::zeroinit_setitem<NumT>::original = funcs->setitem;
		funcs->setitem = &internal::zeroinit_setitem<NumT>::run;
	}

	// Install the guarded numpy dot/inner loop for an MPFR-backed dtype, replacing
	// eigenpy's unguarded SpecialMethods<NumT>::dotfunc.  Call immediately after
	// eigenpy::registerNewType<NumT>().  See internal::guarded_dotfunc.
	template <typename NumT>
	void HardenDotfunc()
	{
		PyArray_Descr *descr = Register::getPyArrayDescr<NumT>();
		PyArray_ArrFuncs *funcs = PyDataType_GetArrFuncs(descr);
		funcs->dotfunc = reinterpret_cast<PyArray_DotFunc*>(&internal::guarded_dotfunc<NumT>);
	}

	// register a single guarded loop on the named numpy ufunc, mirroring the
	// error handling of eigenpy's EIGENPY_REGISTER_*_UFUNC macros.
	inline void registerGuardedLoop(PyObject *numpy, char const *ufunc_name,
	                                int type_code, PyUFuncGenericFunction loop,
	                                int *types, int expected_nargs)
	{
		PyUFuncObject *ufunc =
				(PyUFuncObject *)PyObject_GetAttrString(numpy, ufunc_name);
		if (!ufunc)
		{
			std::stringstream ss;
			ss << "Impossible to define \"" << ufunc_name << "\" for type code "
				 << type_code << std::endl;
			eigenpy::Exception(ss.str());
			return;
		}
		if (expected_nargs != ufunc->nargs)
		{
			PyErr_Format(PyExc_AssertionError,
			             "ufunc %s takes %d arguments, our loop takes %d",
			             ufunc_name, ufunc->nargs, expected_nargs);
			Py_DECREF(ufunc);
			return;
		}
		if (PyUFunc_RegisterLoopForType(ufunc, type_code, loop, types, 0) < 0)
		{
			std::stringstream ss;
			ss << "Impossible to register \"" << ufunc_name << "\" for type code "
				 << type_code << std::endl;
			eigenpy::Exception(ss.str());
		}
		Py_DECREF(ufunc);
	}

	// i lifted this from EigenPy and adapted it: all loops are the guarded
	// versions from internal:: above (eigenpy's read input slots unguarded —
	// see the header comment), and the ordering comparitors are a compile-time
	// option because they are NOT defined for complex types (instantiating
	// them for complex_mp would be a hard error).
	template <typename Scalar, bool WithOrderingComparitors>
	void registerGuardedUfunct()
	{
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
			registerGuardedLoop(numpy, "matmul", type_code,
			                    &internal::guarded_gufunc_matrix_multiply<Scalar>,
			                    types, 3);
		}

		// Binary operators
		{
			int types[3] = {type_code, type_code, type_code};
			registerGuardedLoop(numpy, "add", type_code,
			                    &internal::guarded_binary_op<Scalar, internal::op_add>, types, 3);
			registerGuardedLoop(numpy, "subtract", type_code,
			                    &internal::guarded_binary_op<Scalar, internal::op_subtract>, types, 3);
			registerGuardedLoop(numpy, "multiply", type_code,
			                    &internal::guarded_binary_op<Scalar, internal::op_multiply>, types, 3);
			registerGuardedLoop(numpy, "divide", type_code,
			                    &internal::guarded_binary_op<Scalar, internal::op_divide>, types, 3);
		}

		// Comparison operators
		{
			int types[3] = {type_code, type_code, Register::getTypeCode<bool>()};
			registerGuardedLoop(numpy, "equal", type_code,
			                    &internal::guarded_compare_op<Scalar, internal::op_equal>, types, 3);
			registerGuardedLoop(numpy, "not_equal", type_code,
			                    &internal::guarded_compare_op<Scalar, internal::op_not_equal>, types, 3);

			if constexpr (WithOrderingComparitors) // NOT defined for complex types
			{
				registerGuardedLoop(numpy, "greater", type_code,
				                    &internal::guarded_compare_op<Scalar, internal::op_greater>, types, 3);
				registerGuardedLoop(numpy, "less", type_code,
				                    &internal::guarded_compare_op<Scalar, internal::op_less>, types, 3);
				registerGuardedLoop(numpy, "greater_equal", type_code,
				                    &internal::guarded_compare_op<Scalar, internal::op_greater_equal>, types, 3);
				registerGuardedLoop(numpy, "less_equal", type_code,
				                    &internal::guarded_compare_op<Scalar, internal::op_less_equal>, types, 3);
			}
		}

		// Unary operators
		{
			int types[2] = {type_code, type_code};
			registerGuardedLoop(numpy, "negative", type_code,
			                    &internal::guarded_unary_op<Scalar, internal::op_negative>, types, 2);
			registerGuardedLoop(numpy, "square", type_code,
			                    &internal::guarded_unary_op<Scalar, internal::op_square>, types, 2);
			registerGuardedLoop(numpy, "sqrt", type_code,
			                    &internal::guarded_unary_op<Scalar, internal::op_sqrt>, types, 2);
		}

		Py_DECREF(numpy);
	}

	// kept for call-site compatibility: complex types get no ordering comparitors.
	template <typename Scalar>
	void registerUfunct_without_comparitors()
	{
		registerGuardedUfunct<Scalar, false>();
	}

} // namespace eigenpy

namespace bertini
{
	namespace python
	{

		void EnableEigenPy();

	}
} // namespaces

#endif // include guard
