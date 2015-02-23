

AM_YFLAGS = -d -p `basename $* | sed 's,y$$,,'`
AM_LFLAGS = -s -P`basename $* | sed 's,l$$,,'` -olex.yy.c

AM_CPPFLAGS = -I$(top_srcdir)/include

ACLOCAL_AMFLAGS = -I m4

#set up the empty variables
bin_PROGRAMS =
BUILT_SOURCES =
CLEANFILES =
noinst_LIBRARIES =



include_HEADERS =


SUFFIXES = .c .h .y .l
#see https://www.gnu.org/software/automake/manual/html_node/Suffixes.html



#now include the Makemodule files from the subdirectories.  they will add on to the variables whi

include src/test/classes/Makemodule.am







