

#pragma once

#ifndef BERTINI_PYTHON_EIGENPY_INTERACTION_HPP
#define BERTINI_PYTHON_EIGENPY_INTERACTION_HPP

#include "python_common.hpp"

#include <complex>
#include <cstring>
#include <type_traits>

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


		// template specialization for real numbers.
		//
		// NB: eigenpy's stock getitem (and an earlier version of this one) returns
		// boost::ref(slot) — a Python object ALIASING the numpy buffer.  That is
		// the root of the ADR-0031 hazard family (#259: stored elements silently
		// collapse or SIGABRT once the buffer is reused/freed), and it makes
		// scalars extracted from temporary arrays — np.sum/np.mean results, most
		// visibly — dangle outright.  Returning an owned COPY kills the whole
		// class: an indexed element is a durable value, as numpy users expect.
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
				boost::python::object m(mpfr_scalar); // owned copy — never boost::ref (see above)
				Py_INCREF(m.ptr());
				return m.ptr();
			}
		};

		// a template specialization for complex numbers; see the real one for the
		// copy-not-ref rationale.
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
				boost::python::object m(mpfr_scalar); // owned copy — never boost::ref (see above)
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
		//
		// Writing an output slot: NEVER through BMP operator= on the slot's
		// existing value.  numpy may hand a loop an output slot that bitwise-
		// ALIASES another slot's mpfr allocation -- numpy 2.5 initializes the
		// accumulator of an identityless reduce (np.min/np.max) by memcpy of
		// element 0, so the accumulator and v[0] share one set of limbs.  A
		// plain assignment move-frees the slot's old limbs (freeing v[0]'s
		// storage out from under it: use-after-free, double-free, corrupted
		// allocator, SIGSEGV a few calls later) and writes through shared
		// storage (silently mutating v[0]'s value).  slot_write below is the
		// only sanctioned store: compute the value FIRST (the slot may also be
		// an input), memset the slot to BMP's uninitialized sentinel, then
		// move the fresh value in -- no existing allocation is ever freed or
		// written through.  If the slot held a uniquely-owned value, its
		// allocation leaks (numpy never destructs user-dtype elements anyway);
		// same crash-into-bounded-leak trade as HardenSetitem / ADR-0006.

		// store `val` into a numpy-managed mp slot without freeing or writing
		// through the slot's existing (possibly aliased) allocation.
		template <typename T>
		inline void slot_write(T& slot, T val)
		{
			std::memset(static_cast<void*>(&slot), 0, sizeof(T));
			slot = std::move(val); // move into the sentinel: steals val's limbs, frees nothing
		}

		struct op_add           { template <typename T> static T    apply(T const& x, T const& y) { return T(x + y); } };
		struct op_subtract      { template <typename T> static T    apply(T const& x, T const& y) { return T(x - y); } };
		struct op_multiply      { template <typename T> static T    apply(T const& x, T const& y) { return T(x * y); } };
		struct op_divide        { template <typename T> static T    apply(T const& x, T const& y) { return T(x / y); } };
		// comparison functors are heterogeneous (two type parameters) so the same
		// functor serves the mp-vs-mp loops AND the mp-vs-double tolerance loops
		// (boost::multiprecision compares a number against a double exactly).
		struct op_equal         { template <typename T, typename U> static bool apply(T const& x, U const& y) { return x == y; } };
		struct op_not_equal     { template <typename T, typename U> static bool apply(T const& x, U const& y) { return x != y; } };
		struct op_greater       { template <typename T, typename U> static bool apply(T const& x, U const& y) { return x > y; } };
		struct op_less          { template <typename T, typename U> static bool apply(T const& x, U const& y) { return x < y; } };
		struct op_greater_equal { template <typename T, typename U> static bool apply(T const& x, U const& y) { return x >= y; } };
		struct op_less_equal    { template <typename T, typename U> static bool apply(T const& x, U const& y) { return x <= y; } };

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

		// trait: is this scalar the complex mp type?  several ops (conjugate, sign,
		// the isnan/isinf/isfinite predicates) need a different body for complex.
		template <typename T> struct is_complex_mp : std::false_type {};
		template <> struct is_complex_mp<bertini::complex_mp> : std::true_type {};

		// re-tag `val` to carry `ref`'s precision.  guards against boost's mixed
		// real/complex arithmetic occasionally mis-tagging the result's precision
		// (division is the known offender); .precision(n) preserves the value.
		template <typename T>
		inline T at_precision_of(T val, T const& ref)
		{
			if (val.precision() != ref.precision())
				val.precision(ref.precision());
			return val;
		}

		// ----- unary ops, same-type output ------------------------------------
		// bodies call the same boost::multiprecision free functions the multiprec
		// module binds as scalar functions, so np.exp(arr)[i] == mp.exp(arr[i]).

		struct op_positive   { template <typename T> static T apply(T const& x) { return x; } };
		struct op_reciprocal
		{
			template <typename T> static T apply(T const& x)
			{
				T one(1);
				one.precision(x.precision());
				return at_precision_of(T(one / x), x);
			}
		};
		struct op_conjugate
		{
			template <typename T> static T apply(T const& x)
			{
				if constexpr (is_complex_mp<T>::value)
					return T(conj(x));
				else
					return x;
			}
		};
		// numpy-2 sign semantics: real -> -1/0/+1; complex -> z/|z| (0 at 0).
		struct op_sign
		{
			template <typename T> static T apply(T const& x)
			{
				if constexpr (is_complex_mp<T>::value)
				{
					if (x == 0)
						return at_precision_of(T(0), x);
					bertini::real_mp const mag(abs(x));
					return at_precision_of(T(bertini::real_mp(x.real() / mag),
					                         bertini::real_mp(x.imag() / mag)), x);
				}
				else
				{
					T res(x > 0 ? 1 : (x < 0 ? -1 : 0));
					return at_precision_of(std::move(res), x);
				}
			}
		};

		struct op_exp    { template <typename T> static T apply(T const& x) { return T(exp(x)); } };
		struct op_log    { template <typename T> static T apply(T const& x) { return T(log(x)); } };
		struct op_log10  { template <typename T> static T apply(T const& x) { return T(log10(x)); } };
		struct op_exp2   { template <typename T> static T apply(T const& x) { return T(exp2(x)); } };
		struct op_log2   { template <typename T> static T apply(T const& x) { return T(log2(x)); } };
		struct op_expm1  { template <typename T> static T apply(T const& x) { return T(expm1(x)); } };
		struct op_log1p  { template <typename T> static T apply(T const& x) { return T(log1p(x)); } };
		struct op_cbrt   { template <typename T> static T apply(T const& x) { return T(cbrt(x)); } };

		struct op_sin    { template <typename T> static T apply(T const& x) { return T(sin(x)); } };
		struct op_cos    { template <typename T> static T apply(T const& x) { return T(cos(x)); } };
		struct op_tan    { template <typename T> static T apply(T const& x) { return T(tan(x)); } };
		struct op_arcsin { template <typename T> static T apply(T const& x) { return T(asin(x)); } };
		struct op_arccos { template <typename T> static T apply(T const& x) { return T(acos(x)); } };
		struct op_arctan { template <typename T> static T apply(T const& x) { return T(atan(x)); } };
		struct op_sinh   { template <typename T> static T apply(T const& x) { return T(sinh(x)); } };
		struct op_cosh   { template <typename T> static T apply(T const& x) { return T(cosh(x)); } };
		struct op_tanh   { template <typename T> static T apply(T const& x) { return T(tanh(x)); } };
		struct op_arcsinh { template <typename T> static T apply(T const& x) { return T(asinh(x)); } };
		struct op_arccosh { template <typename T> static T apply(T const& x) { return T(acosh(x)); } };
		struct op_arctanh { template <typename T> static T apply(T const& x) { return T(atanh(x)); } };

		// real-only rounding family
		struct op_floor { template <typename T> static T apply(T const& x) { return T(floor(x)); } };
		struct op_ceil  { template <typename T> static T apply(T const& x) { return T(ceil(x)); } };
		struct op_trunc { template <typename T> static T apply(T const& x) { return T(trunc(x)); } };
		// numpy's rint is round-half-to-EVEN; boost's rint rounds half away from
		// zero, so call mpfr directly in MPFR_RNDN (nearest, ties to even).
		// Unlike the rest of the rounding family, numpy defines rint for complex
		// (component-wise) — np.round on a complex array goes through it.
		struct op_rint
		{
			static bertini::real_mp rint_one(bertini::real_mp const& x)
			{
				bertini::real_mp out(0);
				out.precision(x.precision());
				mpfr_rint(out.backend().data(), x.backend().data(), MPFR_RNDN);
				return out;
			}

			template <typename T> static T apply(T const& x)
			{
				if constexpr (is_complex_mp<T>::value)
					return at_precision_of(T(rint_one(x.real()), rint_one(x.imag())), x);
				else
					return rint_one(x);
			}
		};

		// ----- unary ops, cross-type output -----------------------------------

		// absolute: real -> real, complex -> real (the magnitude).
		struct op_absolute
		{
			template <typename T> static bertini::real_mp apply(T const& x)
			{
				return bertini::real_mp(abs(x));
			}
		};
		struct op_fabs { template <typename T> static T apply(T const& x) { return T(fabs(x)); } };

		// predicates -> bool.  the isnan/isinf/isfinite family are function-like
		// macros in C <math.h>, so call the boost versions qualified.
		struct op_isnan
		{
			template <typename T> static bool apply(T const& x)
			{
				if constexpr (is_complex_mp<T>::value)
					return boost::multiprecision::isnan(x.real()) || boost::multiprecision::isnan(x.imag());
				else
					return boost::multiprecision::isnan(x);
			}
		};
		struct op_isinf
		{
			template <typename T> static bool apply(T const& x)
			{
				if constexpr (is_complex_mp<T>::value)
					return boost::multiprecision::isinf(x.real()) || boost::multiprecision::isinf(x.imag());
				else
					return boost::multiprecision::isinf(x);
			}
		};
		struct op_isfinite
		{
			template <typename T> static bool apply(T const& x)
			{
				if constexpr (is_complex_mp<T>::value)
					return boost::multiprecision::isfinite(x.real()) && boost::multiprecision::isfinite(x.imag());
				else
					return boost::multiprecision::isfinite(x);
			}
		};
		struct op_signbit
		{
			template <typename T> static bool apply(T const& x)
			{
				return boost::multiprecision::signbit(x);
			}
		};

		// ----- binary ops ------------------------------------------------------

		struct op_power    { template <typename T> static T apply(T const& x, T const& y) { return T(pow(x, y)); } };
		struct op_arctan2  { template <typename T> static T apply(T const& x, T const& y) { return T(atan2(x, y)); } };
		struct op_hypot    { template <typename T> static T apply(T const& x, T const& y) { return T(hypot(x, y)); } };
		struct op_copysign { template <typename T> static T apply(T const& x, T const& y) { return T(copysign(x, y)); } };
		// numpy fmod keeps C fmod's sign-of-dividend semantics
		struct op_fmod     { template <typename T> static T apply(T const& x, T const& y) { return T(fmod(x, y)); } };
		// numpy remainder/mod takes the sign of the DIVISOR (python % semantics)
		struct op_remainder
		{
			template <typename T> static T apply(T const& x, T const& y)
			{
				T r(fmod(x, y));
				if (r != 0 && ((r < 0) != (y < 0)))
					r += y;
				return r;
			}
		};
		struct op_floor_divide
		{
			template <typename T> static T apply(T const& x, T const& y)
			{
				return T(floor(x / y));
			}
		};
		// minimum/maximum propagate nan (numpy semantics); fmin/fmax ignore it
		struct op_minimum
		{
			template <typename T> static T apply(T const& x, T const& y)
			{
				if (boost::multiprecision::isnan(x)) return x;
				if (boost::multiprecision::isnan(y)) return y;
				return y < x ? y : x;
			}
		};
		struct op_maximum
		{
			template <typename T> static T apply(T const& x, T const& y)
			{
				if (boost::multiprecision::isnan(x)) return x;
				if (boost::multiprecision::isnan(y)) return y;
				return x < y ? y : x;
			}
		};
		struct op_fmin
		{
			template <typename T> static T apply(T const& x, T const& y)
			{
				if (boost::multiprecision::isnan(x)) return y;
				if (boost::multiprecision::isnan(y)) return x;
				return y < x ? y : x;
			}
		};
		struct op_fmax
		{
			template <typename T> static T apply(T const& x, T const& y)
			{
				if (boost::multiprecision::isnan(x)) return y;
				if (boost::multiprecision::isnan(y)) return x;
				return x < y ? y : x;
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
				slot_write(res, Op::apply(x, y)); // never plain operator= -- see slot_write
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

		// mixed mp-vs-float64 ORDERING comparison (Reversed swaps operand order:
		// false = (mp, double), true = (double, mp)).  Comparisons against a
		// double tolerance -- np.abs(a - b) < 1e-10 -- are safe: boost compares a
		// number against a double exactly, and the result is a bool, so no float
		// ever flows INTO a multiprecision value.  This mirrors the scalar
		// bindings (GreatLessVisitor<T, double>) and the C++ solvers' double
		// ToleranceT.  Deliberately orderings-only: mixed EQUALITY with a float
		// literal is the 0.1-intent trap the unsafe double->mp cast exists to
		// block, and it is not bound at the scalar level either.
		template <typename T, typename Op, bool Reversed>
		void guarded_mixed_compare_op(
				char **args, EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *dimensions,
				EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *steps, void * /*data*/)
		{
			npy_intp is0 = steps[0], is1 = steps[1], os = steps[2], n = *dimensions;
			char *i0 = args[0], *i1 = args[1], *o = args[2];
			const T zero(0);
			for (npy_intp k = 0; k < n; ++k)
			{
				bool& res = *reinterpret_cast<bool*>(o);
				if constexpr (Reversed)
				{
					double const& x = *reinterpret_cast<double const*>(i0);
					T const& y = value_or_zero(*reinterpret_cast<T const*>(i1), zero);
					res = Op::apply(x, y);
				}
				else
				{
					T const& x = value_or_zero(*reinterpret_cast<T const*>(i0), zero);
					double const& y = *reinterpret_cast<double const*>(i1);
					res = Op::apply(x, y);
				}
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
				slot_write(res, Op::apply(x)); // never plain operator= -- see slot_write
				i += is;
				o += os;
			}
		}

		// unary loop with an output type different from the input type
		// (absolute: complex -> real; the isnan family: T -> bool).  Writes into
		// mp-typed output slots go through BMP operator=, which initializes a
		// zeroed destination itself; bool slots are plain bytes.
		template <typename T, typename OutT, typename Op>
		void guarded_unary_op_out(
				char **args, EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *dimensions,
				EIGENPY_NPY_CONST_UFUNC_ARG npy_intp *steps, void * /*data*/)
		{
			npy_intp is = steps[0], os = steps[1], n = *dimensions;
			char *i = args[0], *o = args[1];
			const T zero(0);
			for (npy_intp k = 0; k < n; ++k)
			{
				T const& x = value_or_zero(*reinterpret_cast<T const*>(i), zero);
				OutT& res = *reinterpret_cast<OutT*>(o);
				if constexpr (std::is_trivially_copyable_v<OutT>)
					res = Op::apply(x);
				else
					slot_write(res, Op::apply(x)); // never plain operator= -- see slot_write
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
					slot_write(res, std::move(sum)); // never plain operator= -- see slot_write
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
			slot_write(*reinterpret_cast<T*>(op), std::move(acc)); // never plain operator= -- see slot_write
		}

		// guarded element comparison for the PyArray_ArrFuncs `compare` slot
		// (np.sort / argsort / searchsorted / unique).  eigenpy leaves this slot
		// empty for user dtypes ("type does not have compare function").  Only
		// installed for the real type — complex has no ordering.  nan compares
		// false both ways (weak-ordering violation, same as C doubles): sorting
		// arrays containing nan gives an unspecified nan position, not a crash.
		template <typename T>
		int guarded_compare(const void *a, const void *b, void * /*arr*/)
		{
			const T zero(0);
			T const& x = value_or_zero(*static_cast<T const*>(a), zero);
			T const& y = value_or_zero(*static_cast<T const*>(b), zero);
			if (x < y) return -1;
			if (y < x) return 1;
			return 0;
		}

		// guarded argmax/argmin for the PyArray_ArrFuncs slots (np.argmax /
		// np.argmin / np.max / np.min dispatch through these for user dtypes on
		// some numpy paths).  Mirrors numpy's float semantics: a nan wins
		// immediately (first nan is the arg-extremum).  numpy hands these a
		// contiguous buffer.
		template <typename T, bool Max>
		int guarded_argminmax(void *data, npy_intp n, npy_intp *extremum_ind, void * /*arr*/)
		{
			const T zero(0);
			T const* p = static_cast<T const*>(data);
			*extremum_ind = 0;
			if (n == 0)
				return 0;
			T best = value_or_zero(p[0], zero);
			if (boost::multiprecision::isnan(best))
				return 0;
			for (npy_intp k = 1; k < n; ++k)
			{
				T const& v = value_or_zero(p[k], zero);
				if (boost::multiprecision::isnan(v) || (Max ? best < v : v < best))
				{
					*extremum_ind = k;
					if (boost::multiprecision::isnan(v))
						return 0;
					best = v;
				}
			}
			return 0;
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
			// complex_mp -> real target.  The only such registered cast is
			// complex_mp -> double (complex_mp -> complex128 has its own
			// specialization below; -> integer is unregistered).  Mirror numpy's
			// builtin complex->real cast, which DISCARDS the imaginary part, by
			// converting the real component.  `static_cast<To>(from)` instead
			// routes through boost.multiprecision's complex->scalar conversion,
			// which THROWS std::runtime_error("Could not convert imaginary number
			// to scalar.") whenever the imaginary part is nonzero -- and that C++
			// throw, escaping numpy's C cast loop (an implicitly-noexcept context),
			// calls std::terminate() -> SIGABRT, hard-crashing the interpreter.
			// It bites the moment a complex_mp value with any imaginary part is
			// stored into a real array, e.g. `M = np.zeros(...); M[i,j] = z`.
			return static_cast<To>(from.real());
		}
	};

	// mp -> complex128, for the explicit down-conversion arr.astype(complex)
	// (registered unsafe, like mp -> double: you consciously truncate).
	// boost mp numbers have no conversion operator to std::complex, so go
	// through the components.
	template <>
	struct cast<bertini::real_mp, std::complex<double>>
	{
		static std::complex<double> run(bertini::real_mp const& from)
		{
			if (internal::mpfr_slot<bertini::real_mp>::uninitialized(from))
				return {0.0, 0.0};
			return {from.convert_to<double>(), 0.0};
		}
	};

	template <>
	struct cast<bertini::complex_mp, std::complex<double>>
	{
		static std::complex<double> run(bertini::complex_mp const& from)
		{
			if (internal::mpfr_slot<bertini::complex_mp>::uninitialized(from))
				return {0.0, 0.0};
			return {from.real().convert_to<double>(), from.imag().convert_to<double>()};
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

	// Fill the element-comparison slot (empty in eigenpy's registration), enabling
	// np.sort / np.argsort / np.searchsorted / np.unique.  Real type only —
	// complex has no ordering.  Call immediately after eigenpy::registerNewType.
	template <typename NumT>
	void HardenCompare()
	{
		PyArray_Descr *descr = Register::getPyArrayDescr<NumT>();
		PyArray_ArrFuncs *funcs = PyDataType_GetArrFuncs(descr);
		funcs->compare = reinterpret_cast<PyArray_CompareFunc*>(&internal::guarded_compare<NumT>);
	}

	// Fill the argmax/argmin slots (empty in eigenpy's registration), enabling
	// np.argmax / np.argmin ("data type not ordered" otherwise).  Real type only.
	// Call immediately after eigenpy::registerNewType.
	template <typename NumT>
	void HardenArgMinMax()
	{
		PyArray_Descr *descr = Register::getPyArrayDescr<NumT>();
		PyArray_ArrFuncs *funcs = PyDataType_GetArrFuncs(descr);
		funcs->argmax = reinterpret_cast<PyArray_ArgFunc*>(&internal::guarded_argminmax<NumT, true>);
		funcs->argmin = reinterpret_cast<PyArray_ArgFunc*>(&internal::guarded_argminmax<NumT, false>);
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
	// see the header comment), and the ordering-dependent set is a compile-time
	// option because ordering is NOT defined for complex types (instantiating
	// those functors for complex_mp would be a hard error).  Coverage beyond
	// eigenpy's arithmetic core (absolute, conjugate, the transcendental family,
	// rounding, min/max, the isnan predicates) closes the documented
	// "ufunc not supported" gotchas — every loop body calls the same
	// boost::multiprecision free function the multiprec module binds as the
	// scalar function of the same name.
	template <typename Scalar, bool WithOrderingComparitors>
	void registerGuardedUfunct()
	{
		const int type_code = Register::getTypeCode<Scalar>();
		const int bool_code = Register::getTypeCode<bool>();
		const int real_code = Register::getTypeCode<bertini::real_mp>();

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

		// registration helpers: (in...) -> out signatures.  the types array is
		// copied by PyUFunc_RegisterLoopForType, so stack storage is fine.
		auto unary = [&](char const* name, PyUFuncGenericFunction loop, int out_code)
		{
			int types[2] = {type_code, out_code};
			registerGuardedLoop(numpy, name, type_code, loop, types, 2);
		};
		auto binary = [&](char const* name, PyUFuncGenericFunction loop, int out_code)
		{
			int types[3] = {type_code, type_code, out_code};
			registerGuardedLoop(numpy, name, type_code, loop, types, 3);
		};

		// Matrix multiply
		{
			int types[3] = {type_code, type_code, type_code};
			registerGuardedLoop(numpy, "matmul", type_code,
			                    &internal::guarded_gufunc_matrix_multiply<Scalar>,
			                    types, 3);
		}

		// Binary arithmetic
		binary("add",      &internal::guarded_binary_op<Scalar, internal::op_add>,      type_code);
		binary("subtract", &internal::guarded_binary_op<Scalar, internal::op_subtract>, type_code);
		binary("multiply", &internal::guarded_binary_op<Scalar, internal::op_multiply>, type_code);
		binary("divide",   &internal::guarded_binary_op<Scalar, internal::op_divide>,   type_code);
		binary("power",    &internal::guarded_binary_op<Scalar, internal::op_power>,    type_code);

		// Equality comparisons (defined for real and complex alike)
		binary("equal",     &internal::guarded_compare_op<Scalar, internal::op_equal>,     bool_code);
		binary("not_equal", &internal::guarded_compare_op<Scalar, internal::op_not_equal>, bool_code);

		// Unary, same-type output
		unary("negative",   &internal::guarded_unary_op<Scalar, internal::op_negative>,   type_code);
		unary("positive",   &internal::guarded_unary_op<Scalar, internal::op_positive>,   type_code);
		unary("square",     &internal::guarded_unary_op<Scalar, internal::op_square>,     type_code);
		unary("sqrt",       &internal::guarded_unary_op<Scalar, internal::op_sqrt>,       type_code);
		unary("reciprocal", &internal::guarded_unary_op<Scalar, internal::op_reciprocal>, type_code);
		unary("conjugate",  &internal::guarded_unary_op<Scalar, internal::op_conjugate>,  type_code);
		unary("sign",       &internal::guarded_unary_op<Scalar, internal::op_sign>,       type_code);
		unary("exp",        &internal::guarded_unary_op<Scalar, internal::op_exp>,        type_code);
		unary("log",        &internal::guarded_unary_op<Scalar, internal::op_log>,        type_code);
		unary("log10",      &internal::guarded_unary_op<Scalar, internal::op_log10>,      type_code);
		unary("sin",        &internal::guarded_unary_op<Scalar, internal::op_sin>,        type_code);
		unary("cos",        &internal::guarded_unary_op<Scalar, internal::op_cos>,        type_code);
		unary("tan",        &internal::guarded_unary_op<Scalar, internal::op_tan>,        type_code);
		unary("arcsin",     &internal::guarded_unary_op<Scalar, internal::op_arcsin>,     type_code);
		unary("arccos",     &internal::guarded_unary_op<Scalar, internal::op_arccos>,     type_code);
		unary("arctan",     &internal::guarded_unary_op<Scalar, internal::op_arctan>,     type_code);
		unary("sinh",       &internal::guarded_unary_op<Scalar, internal::op_sinh>,       type_code);
		unary("cosh",       &internal::guarded_unary_op<Scalar, internal::op_cosh>,       type_code);
		unary("tanh",       &internal::guarded_unary_op<Scalar, internal::op_tanh>,       type_code);
		unary("arcsinh",    &internal::guarded_unary_op<Scalar, internal::op_arcsinh>,    type_code);
		unary("arccosh",    &internal::guarded_unary_op<Scalar, internal::op_arccosh>,    type_code);
		unary("arctanh",    &internal::guarded_unary_op<Scalar, internal::op_arctanh>,    type_code);

		// absolute: real -> real, complex -> real (magnitude)
		unary("absolute", &internal::guarded_unary_op_out<Scalar, bertini::real_mp, internal::op_absolute>, real_code);

		// predicates -> bool
		unary("isnan",    &internal::guarded_unary_op_out<Scalar, bool, internal::op_isnan>,    bool_code);
		unary("isinf",    &internal::guarded_unary_op_out<Scalar, bool, internal::op_isinf>,    bool_code);
		unary("isfinite", &internal::guarded_unary_op_out<Scalar, bool, internal::op_isfinite>, bool_code);

		// rint is the one rounding ufunc numpy defines for complex too
		// (component-wise) — np.round dispatches through it
		unary("rint", &internal::guarded_unary_op<Scalar, internal::op_rint>, type_code);

		if constexpr (WithOrderingComparitors) // the ordering-dependent set; NOT defined for complex types
		{
			binary("greater",       &internal::guarded_compare_op<Scalar, internal::op_greater>,       bool_code);
			binary("less",          &internal::guarded_compare_op<Scalar, internal::op_less>,          bool_code);
			binary("greater_equal", &internal::guarded_compare_op<Scalar, internal::op_greater_equal>, bool_code);
			binary("less_equal",    &internal::guarded_compare_op<Scalar, internal::op_less_equal>,    bool_code);

			// mixed mp-vs-float64 orderings, both operand orders: the tolerance
			// idiom `np.abs(a - b) < 1e-10`.  Orderings ONLY -- see
			// guarded_mixed_compare_op for why equality stays mp-vs-mp.
			auto mixed_ordering = [&](char const* name, PyUFuncGenericFunction fwd,
			                          PyUFuncGenericFunction rev)
			{
				int types_td[3] = {type_code, NPY_DOUBLE, bool_code};
				int types_dt[3] = {NPY_DOUBLE, type_code, bool_code};
				registerGuardedLoop(numpy, name, type_code, fwd, types_td, 3);
				registerGuardedLoop(numpy, name, type_code, rev, types_dt, 3);
			};
			mixed_ordering("greater",
			               &internal::guarded_mixed_compare_op<Scalar, internal::op_greater, false>,
			               &internal::guarded_mixed_compare_op<Scalar, internal::op_greater, true>);
			mixed_ordering("less",
			               &internal::guarded_mixed_compare_op<Scalar, internal::op_less, false>,
			               &internal::guarded_mixed_compare_op<Scalar, internal::op_less, true>);
			mixed_ordering("greater_equal",
			               &internal::guarded_mixed_compare_op<Scalar, internal::op_greater_equal, false>,
			               &internal::guarded_mixed_compare_op<Scalar, internal::op_greater_equal, true>);
			mixed_ordering("less_equal",
			               &internal::guarded_mixed_compare_op<Scalar, internal::op_less_equal, false>,
			               &internal::guarded_mixed_compare_op<Scalar, internal::op_less_equal, true>);

			// real-only unary: rounding family (sans rint, registered for both
			// above), fabs, real-only transcendentals
			unary("floor", &internal::guarded_unary_op<Scalar, internal::op_floor>, type_code);
			unary("ceil",  &internal::guarded_unary_op<Scalar, internal::op_ceil>,  type_code);
			unary("trunc", &internal::guarded_unary_op<Scalar, internal::op_trunc>, type_code);
			unary("fabs",  &internal::guarded_unary_op<Scalar, internal::op_fabs>,  type_code);
			unary("exp2",  &internal::guarded_unary_op<Scalar, internal::op_exp2>,  type_code);
			unary("log2",  &internal::guarded_unary_op<Scalar, internal::op_log2>,  type_code);
			unary("expm1", &internal::guarded_unary_op<Scalar, internal::op_expm1>, type_code);
			unary("log1p", &internal::guarded_unary_op<Scalar, internal::op_log1p>, type_code);
			unary("cbrt",  &internal::guarded_unary_op<Scalar, internal::op_cbrt>,  type_code);

			unary("signbit", &internal::guarded_unary_op_out<Scalar, bool, internal::op_signbit>, bool_code);

			// real-only binary
			binary("arctan2",      &internal::guarded_binary_op<Scalar, internal::op_arctan2>,      type_code);
			binary("hypot",        &internal::guarded_binary_op<Scalar, internal::op_hypot>,        type_code);
			binary("copysign",     &internal::guarded_binary_op<Scalar, internal::op_copysign>,     type_code);
			binary("fmod",         &internal::guarded_binary_op<Scalar, internal::op_fmod>,         type_code);
			binary("remainder",    &internal::guarded_binary_op<Scalar, internal::op_remainder>,    type_code);
			binary("floor_divide", &internal::guarded_binary_op<Scalar, internal::op_floor_divide>, type_code);
			binary("minimum",      &internal::guarded_binary_op<Scalar, internal::op_minimum>,      type_code);
			binary("maximum",      &internal::guarded_binary_op<Scalar, internal::op_maximum>,      type_code);
			binary("fmin",         &internal::guarded_binary_op<Scalar, internal::op_fmin>,         type_code);
			binary("fmax",         &internal::guarded_binary_op<Scalar, internal::op_fmax>,         type_code);
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
