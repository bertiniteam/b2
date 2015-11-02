#this is the main Makefile.am for Bertini 2.

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# this project uses non-recursive make.
#
# if you need to add an executable to the core,
# add a Makemodule.am file in the source folder for the main file
# and include it below.
#
# see the b2 github wiki for detailed instructions.
# https://github.com/bertiniteam/b2
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

AM_YFLAGS = -d -p `basename $* | sed 's,y$$,,'`
AM_LFLAGS = -s -P`basename $* | sed 's,l$$,,'` -olex.yy.c

AM_CPPFLAGS = -I$(top_srcdir)/include

ACLOCAL_AMFLAGS = -I m4

####################################
###  set up the empty variables.  these are added onto by the below included files,
###  and are deliberately blank here.
#################################
bin_PROGRAMS =
BUILT_SOURCES =
CLEANFILES =
noinst_LIBRARIES =
EXTRA_PROGRAMS =
include_HEADERS =


#see https://www.gnu.org/software/automake/manual/html_node/Suffixes.html
SUFFIXES = .cpp .hpp



############################
###  now include the Makemodule files from the subdirectories.  they will add on to the variables which are set above
##################################

# include this one first, as it defines some useful groups of files which are used later
# in other files.

include src/bertini/Makemodule.am

include src/classic/Makemodule.am

include src/test/classes/Makemodule.am

include src/test/interactive_testing/Makemodule.am

include src/test/timing/Makemodule.am

include src/test/classic/Makemodule.am

include src/test/tracking_basics/Makemodule.am



