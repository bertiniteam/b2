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

import tracking.amptracking_test as amptracking_test
import tracking.endgame_test as endgame_test
import unittest


mods = (amptracking_test,endgame_test)
suite = unittest.TestSuite();
for tests in mods:
    thissuite = unittest.TestLoader().loadTestsFromModule(tests);
    suite.addTests(thissuite)
#
unittest.TextTestRunner(verbosity=2).run(suite)
