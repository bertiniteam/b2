import zero_dim.basics as basics

import unittest


mods = (basics,)
suite = unittest.TestSuite();
for tests in mods:
    thissuite = unittest.TestLoader().loadTestsFromModule(tests);
    suite.addTests(thissuite)
#
unittest.TextTestRunner(verbosity=2).run(suite)