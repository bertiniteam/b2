#we're building pybertini, version 1.0.alpha4, and the corresponding email is silviana's
AC_INIT([pybertini], [1.0.alpha4], [amethyst@uwec.edu],[pybertini], [http://github.com/bertiniteam/b2])


# Force autoconf to be at least this version number:
AC_PREREQ([2.68])


#
AC_CONFIG_AUX_DIR([config])


# turn on the keeping of produced objects in their folders.  this is for non-recursive make
# and autotools
# see Recursive Make Considered Harmful, and any number of demos.
AM_INIT_AUTOMAKE([1.13 subdir-objects])

#the language for b2 is C++
AC_LANG([C++])


#another directive to use the m4 folder
AC_CONFIG_MACRO_DIR([m4])

#the only produced file will be a single Makefile.
AC_CONFIG_FILES([Makefile])



AX_CXX_COMPILE_STDCXX(17)


#enable the creation of shared libraries
AC_ENABLE_SHARED

# no, do not enable the creation of static libraries
AC_DISABLE_STATIC


#find the linker
AC_PROG_LN_S

#find the command for making directories with their parents
AC_PROG_MKDIR_P

#fire up libtool
LT_INIT

#query python's install location.
AM_PATH_PYTHON([3])


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

AX_PYTHON_DEVEL([>= '3.0.0'])

AX_BOOST_BASE([1.65],, [AC_MSG_ERROR([PyBertini needs Boost at least 1.65, but it was not found in your system])])
AX_BOOST_PYTHON


AX_NUMPY



AX_BOOST_MULTIPRECISION

AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_CHRONO
AX_BOOST_REGEX
AX_BOOST_TIMER
AX_BOOST_UNIT_TEST_FRAMEWORK
AX_BOOST_SERIALIZATION
AX_BOOST_LOG
AX_BOOST_LOG_SETUP
AX_BOOST_THREAD

AX_EIGENPY

AX_BERTINI2



AM_CONFIG_HEADER(config.h)


#wrap it up.
AC_OUTPUT







