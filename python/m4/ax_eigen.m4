#
# SYNOPSIS
#
#   AX_EIGEN()
#
# DESCRIPTION
#
#   Check for the Eigen linear algebra library.
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

AC_DEFUN([AX_EIGEN],
    [
AC_ARG_WITH([eigen],
    [AS_HELP_STRING([--with-eigen@<:@=ARG@:>@], 
        [Use the Eigen linear algebra library from a standard location (ARG=yes),
         from the given location (ARG=<path>),
         or disable it (ARG=no)
         @<:@ARG=yes@:>@ ])],
        [
        if test "$withval" = "no"; then
            want_eigen="no"
        elif test "$withval" = "yes"; then
            want_eigen="yes"
            ac_eigen_path=""
        else
            want_eigen="yes"
            ac_eigen_path="$withval"
        fi
        ],
        [want_eigen="yes"])


if test "x$want_eigen" = "xyes"; then
    AC_SUBST(EIGEN_CPPFLAGS)

    CPPFLAGS_SAVED="$CPPFLAGS"

    AC_REQUIRE([AC_PROG_CXX])
    AC_LANG_PUSH(C++)

    found_eigen_dir=no
    if test "$ac_eigen_path" != ""; then
        
        
        if test -d "$ac_eigen_path/Eigen"; then
            EIGEN_CPPFLAGS="-I$ac_eigen_path"
            found_eigen_dir=yes;
        elif test -d "$ac_eigen_path/eigen3/Eigen"; then
                EIGEN_CPPFLAGS="-I$ac_eigen_path/eigen3/"
                found_eigen_dir=yes;
        elif test -d "$ac_eigen_path/include/eigen3/Eigen"; then
                EIGEN_CPPFLAGS="-I$ac_eigen_path/include/eigen3"
                found_eigen_dir=yes;
        fi

    else 
        for ac_eigen_path_tmp in /usr /usr/local /opt /opt/local /opt/homebrew; do
            if test -d "$ac_eigen_path_tmp/include/eigen3/Eigen"; then
                EIGEN_CPPFLAGS="-I$ac_eigen_path_tmp/include/eigen3"
                found_eigen_dir=yes;
                break;
            elif test -d "$ac_eigen_path_tmp/include/Eigen"; then
                EIGEN_CPPFLAGS="-I$ac_eigen_path_tmp/include"
                found_eigen_dir=yes;
                break;
            fi
        done
    fi

    if test "x$found_eigen_dir" = "xyes"; then
        CPPFLAGS="$CPPFLAGS $EIGEN_CPPFLAGS"
    else
        AC_MSG_ERROR([unable to find Eigen directory])
    fi

    AC_CHECK_HEADERS([Eigen/Dense],
        [
        succeeded=yes;
        CPPFLAGS="$CPPFLAGS_SAVED";
        export CPPFLAGS;
        export EIGEN_CPPFLAGS;
        ],
        [
        CPPFLAGS="$CPPFLAGS_SAVED";
        export CPPFLAGS;
        AC_MSG_ERROR([unable to include Eigen])
        ])
    AC_LANG_POP([C++])
fi

])

       
