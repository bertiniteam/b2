#this is the main Makefile.am for Bertini's Python library.

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# this project uses non-recursive make.
#
# see the b2 github wiki for detailed instructions.
# https://github.com/bertiniteam/b2
# # # # # # # # # # # # # # # # # # # # # # # # # # # #


AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/minieigen/src

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
lib_LTLIBRARIES =
EXTRA_LTLIBRARIES =
lib_LIBRARIES =


#see https://www.gnu.org/software/automake/manual/html_node/Suffixes.html
SUFFIXES = .cpp .hpp



############################
###  now include the Makemodule files from the subdirectories.  they will add on to the variables which are set above
##################################

include src/Makemodule.am

