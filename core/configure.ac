




#we're building b2, version 2.0-alpha6, and the corresponding email is dani brake's
AC_INIT([b2], [2.0-alpha6], [brakeda@uwec.edu],[b2], [http://github.com/bertiniteam/b2])


# Force autoconf to be at least this version number:
AC_PREREQ([2.68])


#
AC_CONFIG_AUX_DIR([config])

# turn on the keeping of produced objects in their folders.  this is for non-recursive make
# and autotools
# see Recursive Make Considered Harmful, and any number of demos.
AM_INIT_AUTOMAKE([1.13 subdir-objects])

#the language for Bertini2 is C++
AC_LANG([C++])


#another directive to use the m4 folder
AC_CONFIG_MACRO_DIR([m4])

#the only produced file will be a single Makefile.
AC_CONFIG_FILES([Makefile])


#find the flags for C++14
AX_CXX_COMPILE_STDCXX_14


#enable the creation of shared libraries
AC_ENABLE_SHARED

#enable the creation of static libraries
AC_ENABLE_STATIC

#find the linker
AC_PROG_LN_S

#find the command for making directories with their parents
AC_PROG_MKDIR_P

#fire up libtool
LT_INIT



#AC_ARG_ENABLE(option-name, help-string, action-if-present, action-if-not-present)


AC_ARG_ENABLE([forbid_mixed_arithmetic],
    AS_HELP_STRING([--disable-forbid_mixed_arithmetic], [Disable the forbidding of mixed type arithmetic between built-in floating point types, and fixed or variable multiple precision types (mpfr, etc).  By default, it is forbidden at compile time to mix doubles and multiprecisions.  If you want to be able to multiply a 60-digit mpfr_float by a 16-digit double (and introduce error around the 17th digit of the double), then disable this.  Feels like a double-negative to disable forbiddence, but hey, it's more unsafe to mix doubles and mpfr_floats than it is to speak a double negative.  Have a nice day!]))

AS_IF([test "x$enable_forbid_mixed_arithmetic" != "xno"],[
  AC_DEFINE([B2_FORBID_MIXED_ARITHMETIC], [1],[Define if B2_FORBID_MIXED_ARITHMETIC keyword should be used to not forbid mixed double-mpfr arithmetic.  ])
])


AC_ARG_ENABLE([thread_local],
    AS_HELP_STRING([--disable-thread_local], [Disable thread_local storage for fast complex multiple-precision arithmetic, and other places.  XCode coming with OSX Mavericks and before do NOT implement this keyword.  Since this keyword provides the thread-safety mechanism for complex arithmetic, disabling thread_local threatens thread safety.  Disable at your own risk.  Moving to a newer OSX and XCode to get `thread_local` is the best solution.]))

AS_IF([test "x$enable_thread_local" != "xno"],[
	AC_DEFINE([USE_THREAD_LOCAL], [1],[Define if thread_local keyword should be used.  ])
])


AC_ARG_ENABLE([expression_templates],
    AS_HELP_STRING([--disable-expression_templates], [Disable the use of expression templates, particularly those coming from Boost.Multiprecision.]))

AS_IF([test "x$enable_expression_templates" != "xno"],[
	AC_DEFINE([BMP_EXPRESSION_TEMPLATES], [1],[Use expression templates from Boost.Multiprecision.])
])


AC_ARG_ENABLE([copy_and_swap],
    AS_HELP_STRING([--enable-copy_and_swap], [Enable the use of copy-and-swap for exception safety of the library.]),
    [],
    [enable_copy_and_swap=no])

AS_IF([test "x$enable_copy_and_swap" != "xno"],[
	AC_DEFINE([BERTINI_ENABLE_COPY_AND_SWAP], [1],[Enable copy-and-swap idiom in fundamental types, for exception safety.])
])


AC_ARG_ENABLE([unity_build],
[  --enable-unity_build    Turn on unity_build-style building, better for a low-thread environment.  possible values yes and no.],
[case "${enableval}" in
  yes) unity_build=true ;;
  no)  unity_build=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-unity_build]) ;;
esac],[unity_build=false])
AM_CONDITIONAL([UNITYBUILD], [test x$unity_build = xtrue])


# the form of the following commands --
# AC_SEARCH_LIBS(function, libraries-list, action-if-found, action-if-not-found, extra-libraries)


#find the math library
AC_SEARCH_LIBS([cos], [m], [], [
  AC_MSG_ERROR([unable to find the cos() function])
  ])

#find gmp
AC_SEARCH_LIBS([__gmpz_init],[gmp], [],[
	AC_MSG_ERROR([unable to find gmp])
])

# find mpfr
AC_SEARCH_LIBS([mpfr_get_version],[mpfr], [],[
  AC_MSG_ERROR([unable to find mpfr])
])



# look for a header file in Eigen, and croak if fail to find.
AX_EIGEN

AX_BOOST_BASE([1.62],, [AC_MSG_ERROR([Bertini2 needs Boost at least 1.62 but it was not found in your system])])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_CHRONO
AX_BOOST_REGEX
AX_BOOST_TIMER
AX_BOOST_UNIT_TEST_FRAMEWORK
AX_BOOST_SERIALIZATION
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_LOG
AX_BOOST_LOG_SETUP
AX_BOOST_THREAD




AM_CONFIG_HEADER(include/bertini2/config.h)


#wrap it up.
AC_OUTPUT
