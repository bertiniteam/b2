#
# SYNOPSIS
#
#   AX_BERTINI2()
#
# DESCRIPTION
#
#   Check for the Bertini2 numerical algebraic geometry library.
#
#
# LICENSE
#
#   Copyright Daniel Brake 2016
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

AC_DEFUN([AX_BERTINI2],
    [
AC_ARG_WITH([bertini],
    [AS_HELP_STRING([--with-bertini@<:@=ARG@:>@], 
        [Use the Bertini2 numerical algebraic geometry library from a specified or standard location (ARG=yes),
         from the given location (ARG=<path>),
         or disable it (ARG=no)
         @<:@ARG=yes@:>@ ])],
        [
        if test "$withval" = "no"; then
            want_bertini="no"
        elif test "$withval" = "yes"; then
            want_bertini="yes"
            ac_bertini_path=""
        else
            want_bertini="yes"
            ac_bertini_path="$withval"
        fi
        ],
        [want_bertini="yes"])


if test "x$want_bertini" = "xyes"; then

    AC_REQUIRE([AC_PROG_CXX])
    AC_LANG_PUSH(C++)

    CPPFLAGS_SAVED="$CPPFLAGS"
    found_bertini_include_dir=no
    if test "$ac_bertini_path" != ""; then
        
        if test -d "$ac_bertini_path/bertini2"; then
            BERTINI_CPPFLAGS="-I$ac_bertini_path"
            found_bertini_include_dir=yes;
        elif test -d "$ac_bertini_path/include/bertini2"; then
            BERTINI_CPPFLAGS="-I$ac_bertini_path/include"
            found_bertini_include_dir=yes;
        fi

    else 
        for ac_bertini_path_tmp in /usr /usr/local /opt /opt/local ; do
            if test -d "$ac_bertini_path_tmp/include/bertini2"; then
                BERTINI_CPPFLAGS="-I$ac_bertini_path_tmp/include"
                found_bertini_include_dir=yes;
                break;
            fi
        done
    fi

    if test "x$found_bertini_include_dir" = "xyes"; then
        CPPFLAGS="$CPPFLAGS_SAVED $BERTINI_CPPFLAGS"
        export CPPFLAGS
        AC_MSG_NOTICE([found bertini2 includedir])
    else
        AC_MSG_ERROR([unable to find include/bertini2 directory])
    fi

    LDFLAGS_SAVED="$LDFLAGS"
    found_bertini_lib_dir=no
    if test "$ac_bertini_path" != ""; then
        if test -d "$ac_bertini_path/lib"; then
            BERTINI_LDFLAGS="-L$ac_bertini_path/lib"
            found_bertini_lib_dir=yes;
        fi
    else 
        for ac_bertini_path_tmp in /usr /usr/local /opt /opt/local ; do
            if test -d "$ac_bertini_path_tmp/lib"; then
                BERTINI_LDFLAGS="-L$ac_bertini_path_tmp/lib"
                found_bertini_lib_dir=yes;
                break;
            fi
        done
    fi

    if test "x$found_bertini_lib_dir" = "xyes"; then
        LDFLAGS="$LDFLAGS_SAVED $BERTINI_LDFLAGS"
        export LDFLAGS
        AC_MSG_NOTICE([found bertini2 libdir])
    else
        AC_MSG_ERROR([unable to find lib directory for bertini2])
    fi

    AC_CHECK_HEADERS([bertini2/mpfr_complex.hpp],
        [
            succeeded=yes;
            AC_SUBST(BERTINI_CPPFLAGS)
        ],
        [
        AC_MSG_ERROR([unable to find required file mpfr_complex.hpp in include for Bertini2])
        ])

    AC_SEARCH_LIBS(
        [HaveBertini2],
        [bertini2], 
        [
            AC_SUBST(BERTINI_LDFLAGS)
        ],
        [
        AC_MSG_ERROR([unable to find bertini2's library])],
        [$BOOST_LDFLAGS]
       )

    AC_LANG_POP([C++])

    CPPFLAGS="$CPPFLAGS_SAVED";
    LDFLAGS="$LDFLAGS_SAVED";
fi

])

       
