# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_boost_timer.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_TIMER
#
# DESCRIPTION
#
#   Test for Timer library from the Boost C++ libraries. The macro requires
#   a preceding call to AX_BOOST_BASE. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_TIMER_LIB)
#
#   And sets:
#
#     HAVE_BOOST_TIMER
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Michael Tindal
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 22

AC_DEFUN([AX_BOOST_TIMER],
[
AC_ARG_WITH([boost-timer],
AS_HELP_STRING([--with-boost-timer@<:@=special-lib@:>@],
[use the Timer library from boost - it is possible to specify a certain library for the linker
e.g. --with-boost-timer=boost_timer-gcc-mt-d-1_33_1 ]),
[
if test "$withval" = "no"; then
want_boost="no"
elif test "$withval" = "yes"; then
want_boost="yes"
ax_boost_user_timer_lib=""
else
want_boost="yes"
ax_boost_user_timer_lib="$withval"
fi
],
[want_boost="yes"]
)

if test "x$want_boost" = "xyes"; then
AC_REQUIRE([AC_PROG_CC])
CPPFLAGS_SAVED="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
export CPPFLAGS

LDFLAGS_SAVED="$LDFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
export LDFLAGS

AC_CACHE_CHECK(whether the Boost::Timer library is available,
ax_cv_boost_timer,
[AC_LANG_PUSH([C++])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/timer.hpp>
]],
[[boost::timer r(); return 0;]])],
ax_cv_boost_timer=yes, ax_cv_boost_timer=no)
AC_LANG_POP([C++])
])
if test "x$ax_cv_boost_timer" = "xyes"; then
AC_DEFINE(HAVE_BOOST_TIMER,,[define if the Boost::Timer library is available])
BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
if test "x$ax_boost_user_timer_lib" = "x"; then
for libextension in `ls $BOOSTLIBDIR/libboost_timer*.so* $BOOSTLIBDIR/libboost_timer*.dylib* $BOOSTLIBDIR/libboost_timer*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_timer.*\)\.so.*$;\1;' -e 's;^lib\(boost_timer.*\)\.dylib.*;\1;' -e 's;^lib\(boost_timer.*\)\.a.*$;\1;'` ; do
ax_lib=${libextension}
AC_CHECK_LIB($ax_lib, exit,
[BOOST_TIMER_LIB="-l$ax_lib"; AC_SUBST(BOOST_TIMER_LIB) link_timer="yes"; break],
[link_timer="no"])
done
if test "x$link_timer" != "xyes"; then
for libextension in `ls $BOOSTLIBDIR/boost_timer*.dll* $BOOSTLIBDIR/boost_timer*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_timer.*\)\.dll.*$;\1;' -e 's;^\(boost_timer.*\)\.a.*$;\1;'` ; do
ax_lib=${libextension}
AC_CHECK_LIB($ax_lib, exit,
[BOOST_TIMER_LIB="-l$ax_lib"; AC_SUBST(BOOST_TIMER_LIB) link_timer="yes"; break],
[link_timer="no"])
done
fi

else
for ax_lib in $ax_boost_user_timer_lib boost_timer-$ax_boost_user_timer_lib; do
AC_CHECK_LIB($ax_lib, main,
[BOOST_TIMER_LIB="-l$ax_lib"; AC_SUBST(BOOST_TIMER_LIB) link_timer="yes"; break],
[link_timer="no"])
done
fi
if test "x$ax_lib" = "x"; then
AC_MSG_ERROR(Could not find a version of the Boost::Timer library!)
fi
if test "x$link_timer" != "xyes"; then
AC_MSG_ERROR(Could not link against $ax_lib !)
fi
fi

CPPFLAGS="$CPPFLAGS_SAVED"
LDFLAGS="$LDFLAGS_SAVED"
fi
])
