#
# SYNOPSIS
#
#   AX_NUMPY()
#
# DESCRIPTION
#
#   Check for the Numpy linear algebra library.
#
#
# LICENSE
#
#   Copyright Silviana Amethyst 2016-2018
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_NUMPY],
    [
AC_ARG_WITH([numpy],
    [AS_HELP_STRING([--with-numpy@<:@=ARG@:>@], 
        [Use numpy from the automatically determined location (ARG=yes),
         from the given location (ARG=<path>),
         or disable it (ARG=no)
         @<:@ARG=yes@:>@ ])],
        [
        if test "$withval" = "no"; then
            want_numpy="no"
        elif test "$withval" = "yes"; then
            want_numpy="yes"
            ac_numpy_path=""
        else
            want_numpy="yes"
            ac_numpy_path="$withval"
        fi
        ],
        [want_numpy="yes"])




if test "x$want_numpy" = "xyes"; then


  CPPFLAGS_SAVED="$CPPFLAGS"
  AC_REQUIRE([AC_PROG_CXX])
  AC_LANG_PUSH(C++)



  AC_REQUIRE([AM_PATH_PYTHON])
  AC_MSG_CHECKING([for numpy library can import])
  $PYTHON -c "import numpy" 2>/dev/null
  if test $? == 0; then
    AC_MSG_RESULT([yay, able to `import numpy`])
  else
    AC_MSG_FAILURE([unable to call `import numpy` in your Python install.  I conclude: numpy library not found for your Python installation.  please install numpy via your package manager])
  fi



  found_numpy_dir=no
  if test "$ac_numpy_path" != ""; then
      AC_MSG_CHECKING([])
      if test -d "$ac_numpy_path/numpy/"; then
          NUMPY_CPPFLAGS="-I$ac_numpy_path"
          found_numpy_dir=yes;
      fi

  else 
      AC_CACHE_CHECK([for numpy include directory],
          [_b2_numpy_incdir],
          [
            _b2_numpy_incdir=`$PYTHON -c "import numpy as np; numpy_include_path=np.get_include(); print(numpy_include_path)"`; found_numpy_dir=yes;
          ])
      AC_SUBST([NUMPY_CPPFLAGS], [-I$_b2_numpy_incdir])

  fi


    if test "x$found_numpy_dir" = "xyes"; then
        CPPFLAGS="$CPPFLAGS $NUMPY_CPPFLAGS"
    else
        AC_MSG_ERROR([unable to find headers for numpy.])
    fi

    AC_CHECK_HEADERS([numpy/numpyconfig.h],
        [
        succeeded=yes;
        CPPFLAGS="$CPPFLAGS_SAVED";
        export CPPFLAGS;
        export NUMPY_CPPFLAGS;
        ],
        [
        CPPFLAGS="$CPPFLAGS_SAVED";
        export CPPFLAGS;
        AC_MSG_ERROR([unable to include numpy.  see config.log for details])
        ])
    AC_LANG_POP([C++])
fi

])


