#
# SYNOPSIS
#
#   AX_BOOST_MULTIPRECISION()
#
# DESCRIPTION
#
#   Check for the Boost Multiprecision library, allowing for it to be separate from the install of Boost.
#
#
# LICENSE
#
#   Copyright Silviana Amethyst, 2021
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_BOOST_MULTIPRECISION],
    [
AC_ARG_WITH([boost_multiprecision],
    [AS_HELP_STRING([--with-boost_multiprecision@<:@=ARG@:>@], 
        [Use the Boost multiprecision library from the standard location (ARG=yes),
         from a specific location (ARG=<path>)
         @<:@ARG=yes@:>@ ])],
        [
        if test "$withval" = "no"; then
            want_boost_multiprecision="no"
        elif test "$withval" = "yes"; then
            want_boost_multiprecision="yes"
            ac_boost_multiprecision_path=""
        else
            want_boost_multiprecision="yes"
            ac_boost_multiprecision_path="$withval"
        fi
        ],
        [want_boost_multiprecision="yes"])


if test "x$want_boost_multiprecision" = "xyes"; then
    
    CPPFLAGS_SAVED="$CPPFLAGS"

    AC_REQUIRE([AC_PROG_CXX])
    AC_LANG_PUSH(C++)

    found_bmp_dir=no
    if test "$ac_boost_multiprecision_path" != ""; then

        
        AC_MSG_CHECKING([for Boost.Multiprecision at the path you specified: $ac_boost_multiprecision_path, by looking for $ac_boost_multiprecision_path/boost/multiprecision/mpc.hpp])
        if test -f "$ac_boost_multiprecision_path/boost/multiprecision/mpc.hpp"; then
            BOOST_MULTIPRECISION_CPPFLAGS="-I$ac_boost_multiprecision_path";
            found_bmp_dir=yes;
        else
            AC_MSG_CHECKING([for Boost.Multiprecision at the path you specified, by looking for $ac_boost_multiprecision_path/include/boost/multiprecision/mpc.hpp])
            if test -f "$ac_boost_multiprecision_path/include/boost/multiprecision/mpc.hpp"; then
                BOOST_MULTIPRECISION_CPPFLAGS="-I$ac_boost_multiprecision_path/include";
                found_bmp_dir=yes;
            fi
        fi

    else
        AC_MSG_CHECKING([for Boost.Multiprecision as part of the found Boost installation: $BOOST_CPPFLAGS])

        CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
        export CPPFLAGS
        AC_COMPILE_IFELSE([
             AC_LANG_PROGRAM(
                 [[@%:@include <boost/multiprecision/cpp_dec_float.hpp>]],
                 [[using boost::multiprecision::cpp_dec_float_50;
                   return 0;]])],
             found_bmp_dir=yes, found_bmp_dir=no)

        CPPFLAGS="$CPPFLAGS_SAVED"
        export CPPFLAGS

    fi

    if test "x$found_bmp_dir" = "xyes"; then
        BOOST_CPPFLAGS="$BOOST_MULTIPRECISION_CPPFLAGS $BOOST_CPPFLAGS"
    else
        AC_MSG_ERROR([unable to find Boost Multiprecision.  See `config.log`.])
    fi

    CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
    export CPPFLAGS
    AC_CHECK_HEADERS([boost/multiprecision/mpc.hpp],
        [
        succeeded=yes;
        CPPFLAGS="$CPPFLAGS_SAVED";
        export CPPFLAGS;
        ],
        [
        CPPFLAGS="$CPPFLAGS_SAVED";
        export CPPFLAGS;
        AC_MSG_ERROR([unable to include Boost.Multiprecision])
        ])
    AC_LANG_POP([C++])
else
    AC_MSG_ERROR([you said you don't want Boost.Multiprecision, but it's required for Bertini 2.])
fi

])

      