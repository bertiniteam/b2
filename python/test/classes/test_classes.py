# This file is part of Bertini 2.
# 
# python/test/b2_class_test.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# python/test/b2_class_test.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with python/test/b2_class_test.py.  If not, see <http://www.gnu.org/licenses/>.
# 
#  Copyright(C) 2016 by Bertini2 Development Team
# 
#  See <http://www.gnu.org/licenses/> for a copy of the license, 
#  as well as COPYING.  Bertini2 is provided with permitted 
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
# 
#   James Collins
#   West Texas A&M University
#   Spring 2016
# 

from __future__ import print_function

import classes.mpfr_test as mpfr_test
import classes.function_tree_test as function_tree_test
import classes.differentiation_test as differentiation_test
import classes.system_test as system_test
import classes.parser_test as parser_test

import unittest




mods = (mpfr_test, function_tree_test, differentiation_test, system_test, parser_test)
suite = unittest.TestSuite();
print(mods)
for tests in mods:
    thissuite = unittest.TestLoader().loadTestsFromModule(tests);
    print(thissuite)
    suite.addTests(thissuite)
#
unittest.TextTestRunner(verbosity=2).run(suite)
