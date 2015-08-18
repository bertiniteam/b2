




#we're building b2, version 2.0.0, and the corresponding email is dan brake's
AC_INIT([b2], [2015.06.01.dev], [danielthebrake@gmail.com],[b2], [http://github.com/bertiniteam/b2])


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





#this calls a file in the m4/ directory, which sets up the MPI wrapper stuffs
LX_FIND_MPI

# print a notice about the presence of MPI on the machine.
AS_IF([test "$have_CXX_mpi" = "yes"],
	AC_MSG_NOTICE([mpi c++ appears to work correctly]),
	AC_MSG_NOTICE([it appears you do not have mpi installed])
	)



# the form of the following commands --
# AC_SEARCH_LIBS(function, libraries-list, action-if-found, action-if-not-found, extra-libraries)

# find mpfr
AC_SEARCH_LIBS([mpfr_get_version],[mpfr], [],[
	AC_MSG_ERROR([unable to find mpfr])
])

#find gmp
AC_SEARCH_LIBS([__gmpz_init],[gmp], [],[
	AC_MSG_ERROR([unable to find gmp])
])


#find the math library
AC_SEARCH_LIBS([cos], [m], [], [
  AC_MSG_ERROR([unable to find the cos() function])
  ])

# look for a header file in Eigen, and croak if fail to find.
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([eigen3/Eigen/Dense],[],AC_MSG_ERROR([unable to find Eigen 3  (eigen3/Eigen/Dense)]),[])
AC_LANG_POP([C++])


AX_BOOST_BASE([1.55],, [AC_MSG_ERROR([bertini2 needs Boost at least 1.55, but it was not found in your system])])
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_CHRONO
AX_BOOST_REGEX
AX_BOOST_TIMER
AX_BOOST_UNIT_TEST_FRAMEWORK
AX_BOOST_SERIALIZATION






AM_CONFIG_HEADER(config.h)


#wrap it up.
AC_OUTPUT







