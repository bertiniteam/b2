#
# SYNOPSIS
#
#   AX_EIGENPY()
#
# DESCRIPTION
#
#   Find the library for EigenPy
#
#
# LICENSE
#
#   Copyright silviana amethyst 2023
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_EIGENPY],
    [
AC_ARG_WITH([eigenpy],
    [AS_HELP_STRING([--with-eigenpy@<:@=ARG@:>@],
        [Use eigenpy from a specified or standard location (ARG=yes),
         from the given location (ARG=<path>),
         or disable it (ARG=no)
         @<:@ARG=yes@:>@ ])],
        [
        if test "$withval" = "no"; then
            want_eigenpy="no"
        elif test "$withval" = "yes"; then
            want_eigenpy="yes"
            ac_eigenpy_path=""
        else
            want_eigenpy="yes"
            ac_eigenpy_path="$withval"
        fi
        ],
        [want_eigenpy="yes"])


if test "x$want_eigenpy" = "xyes"; then

    AC_REQUIRE([AC_PROG_CXX])
    AC_LANG_PUSH(C++)

    CPPFLAGS_SAVED="$CPPFLAGS"
    LDFLAGS_SAVED="$LDFLAGS"

    found_eigenpy_include_dir=no
    found_eigenpy_lib_dir=no
    if test "$ac_eigenpy_path" != ""; then

        if test -d "$ac_eigenpy_path/include/eigenpy"; then

            if test -d "$ac_eigenpy_path/lib"; then
                EIGENPY_CPPFLAGS="-I$ac_eigenpy_path/include";
                EIGENPY_LDFLAGS="-L$ac_eigenpy_path/lib -leigenpy";
                found_eigenpy_include_dir=yes;
                found_eigenpy_lib_dir=yes;
            else
                AC_MSG_ERROR([unable to find lib directory for eigenpy symmetric to include directory $ac_eigenpy_path/include/eigenpy])
            fi
        else
            AC_MSG_ERROR([unable to find include/eigenpy directory at specified location --with-eigenpy=$ac_eigenpy_path])
        fi

    else
        for ac_eigenpy_path_tmp in /usr /usr/local /opt /opt/local /opt/homebrew ; do
            if test -d "$ac_eigenpy_path_tmp/include/eigenpy"; then
                if test -d "$ac_eigenpy_path_tmp/lib"; then
                    found_eigenpy_include_dir=yes;
                    found_eigenpy_lib_dir=yes;
                    EIGENPY_LDFLAGS="-L$ac_eigenpy_path_tmp/lib -leigenpy";
                    EIGENPY_CPPFLAGS="-I$ac_eigenpy_path_tmp/include";
                    break;
                else
                    AC_MSG_ERROR([unable to find lib directory for eigenpy symmetric to include directory $ac_eigenpy_path_tmp/include/eigenpy"])
                fi
            fi
        done
        if test "x$found_eigenpy_include_dir" = "xno"; then
            AC_MSG_ERROR([unable to find include/eigenpy directory at standard locations (/usr /usr/local /opt /opt/local /opt/homebrew)])
        fi
    fi

    found_eigenpy_dir=no
    if test "x$found_eigenpy_include_dir" = "xyes"; then
      if test "x$found_eigenpy_lib_dir" = "xyes"; then
          found_eigenpy_dir=yes;
      fi
    fi


    if test "x$found_eigenpy_dir" = "xyes"; then
        CPPFLAGS="$CPPFLAGS_SAVED $BOOST_CPPFLAGS $EIGENPY_CPPFLAGS";
        LDFLAGS="$LDFLAGS_SAVED $BOOST_LDFLAGS $EIGENPY_LDFLAGS";
        export CPPFLAGS;
        export LDFLAGS;
    else
        AC_MSG_ERROR([unable to find eigenpy])
    fi


    AC_CHECK_HEADERS([eigenpy/version.hpp],
        [
            succeeded=yes;
            AC_SUBST(EIGENPY_CPPFLAGS)
        ],
        [
        AC_MSG_ERROR([unable to find required file version.hpp in include for eigenpy])
        ])

    AC_SUBST(EIGENPY_LDFLAGS)


    AC_MSG_NOTICE([found eigenpy.  EIGENPY_CPPFLAGS=$EIGENPY_CPPFLAGS   EIGENPY_LDFLAGS=$EIGENPY_LDFLAGS])

    AC_LANG_POP([C++])

    CPPFLAGS="$CPPFLAGS_SAVED";
    LDFLAGS="$LDFLAGS_SAVED";
fi

])
