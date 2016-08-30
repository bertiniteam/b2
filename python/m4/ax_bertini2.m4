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
AC_ARG_WITH([bertini2],
    [AS_HELP_STRING([--with-bertini2@<:@=ARG@:>@], 
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
    
    CPPFLAGS_SAVED="$CPPFLAGS"

    AC_REQUIRE([AC_PROG_CXX])
    AC_LANG_PUSH(C++)

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
        CPPFLAGS="$CPPFLAGS $BERTINI_CPPFLAGS"
    else
        AC_MSG_ERROR([unable to find include/bertini2 directory])
    fi


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
        LDFLAGS="$CPPFLAGS $BERTINI_LDFLAGS"
    else
        AC_MSG_ERROR([unable to find lib directory for bertini2])
    fi

    AC_CHECK_HEADERS([bertini2/mpfr_complex.hpp],
        [
        succeeded=yes;
        export CPPFLAGS;
        ],
        [
        CPPFLAGS="$CPPFLAGS_SAVED";
        export CPPFLAGS;
        AC_MSG_ERROR([unable to find required file mpfr_complex.hpp in include for Bertini2])
        ])

    AC_SEARCH_LIBS([HaveBertini2],[bertini2], [],
       [
        AC_MSG_ERROR([unable to find bertini2's library])],
       [$BOOST_FILESYSTEM_LIB $BOOST_SYSTEM_LIB $BOOST_CHRONO_LIB $BOOST_REGEX_LIB $BOOST_TIMER_LIB $BOOST_UNIT_TEST_FRAMEWORK_LIB $BOOST_SERIALIZATION_LIB $MPI_CXXLDFLAGS ])

    AC_LANG_POP([C++])
fi

])

       
