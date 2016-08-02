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

EXTRA_PROGRAMS =
EXTRA_LTLIBRARIES =

include_HEADERS =

lib_LTLIBRARIES =
lib_LIBRARIES =

noinst_LIBRARIES =
noinst_LTLIBRARIES =

TESTS =

######  end deliberately blank items

#see https://www.gnu.org/software/automake/manual/html_node/Suffixes.html
SUFFIXES = .cpp .hpp

rootincludedir = $(includedir)/bertini2
rootinclude_HEADERS = include/bertini2/config.h

#initialize to empty and add to it in the Makemodule.am files below


############################
###  now include the Makemodule files from the subdirectories.  they will add on to the variables which are set above
##################################

# include the bertini Makemodule first, as it defines some useful groups of files which are used later
# in other files.

include src/bertini/Makemodule.am
include src/basics/Makemodule.am
include src/function_tree/Makemodule.am
include src/system/Makemodule.am
include src/tracking/Makemodule.am
include src/detail/Makemodule.am

include test/classes/Makemodule.am
include test/tracking_basics/Makemodule.am
include test/classic/Makemodule.am
include test/settings/Makemodule.am
include test/timing/Makemodule.am
include test/endgames/Makemodule.am
include test/pools/Makemodule.am
include test/generating/Makemodule.am
include test/nag_algorithms/Makemodule.am

