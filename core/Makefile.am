#this is the main Makefile.am for Bertini 2.

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# this project uses non-recursive make.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

AM_YFLAGS = -d -p `basename $* | sed 's,y$$,,'`
AM_LFLAGS = -s -P`basename $* | sed 's,l$$,,'` -olex.yy.c

AM_CPPFLAGS = -I$(top_srcdir)/include $(BOOST_CPPFLAGS)

ACLOCAL_AMFLAGS = -I m4

####################################
###  set up the empty variables.  these are added onto by the below included files,
###  and are deliberately blank here.
#################################

bin_PROGRAMS =
BUILT_SOURCES =
CLEANFILES =

#installed headers
include_HEADERS =
rootinclude_HEADERS =

#installed libraries, both free-standing and libtool
lib_LTLIBRARIES =
lib_LIBRARIES =

#helper libraries, both free-standing and libtool
noinst_LIBRARIES =
noinst_LTLIBRARIES =

#programs which are not installed, but are optional targets for building.
EXTRA_PROGRAMS =
EXTRA_LTLIBRARIES =

TESTS =

core_all =
core_sources =
core_headers =


############################################
######  end deliberately blank items     ##
##########################################


#see https://www.gnu.org/software/automake/manual/html_node/Suffixes.html
SUFFIXES = .cpp .hpp

rootincludedir = $(includedir)/bertini2

#initialize to empty and add to it in the Makemodule.am files below

############
#
#  a note for developers:
#
######
#
# if you need to add an executable to the core,
# add a Makemodule.am file in the source folder for the `main()` file
# and include it below.
#
# see the b2 github wiki for detailed instructions.
# https://github.com/bertiniteam/b2
##############

############################
###  now include the Makemodule.am files from the subdirectories;
###  they will add on to the variables which are set above.
###  note that the name `Makemodule.am` is arbitrary.
##################################

# include the bertini2 Makemodule first, as it defines some useful groups of files which are used later
# in other files.

include src/basics/Makemodule.am
include src/common/Makemodule.am
include src/function_tree/Makemodule.am
include src/system/Makemodule.am
include src/tracking/Makemodule.am
include src/endgames/Makemodule.am
include src/detail/Makemodule.am
include src/io/Makemodule.am
include src/nag_algorithms/Makemodule.am
include src/nag_datatypes/Makemodule.am
include src/pool/Makemodule.am
include src/parallel/Makemodule.am

include src/corelibrary/Makemodule.am
include src/blackbox/Makemodule.am

###
#  and finally test suites, built as extras
####

include test/classes/Makemodule.am
include test/tracking_basics/Makemodule.am
include test/classic/Makemodule.am
include test/settings/Makemodule.am
include test/endgames/Makemodule.am
include test/pools/Makemodule.am
include test/generating/Makemodule.am
include test/nag_algorithms/Makemodule.am
include test/nag_datatypes/Makemodule.am
include test/blackbox/Makemodule.am
