import nid.nid_test as nid_test

import unittest


mods = (nid_test,)
suite = unittest.TestSuite();
for tests in mods:
    thissuite = unittest.TestLoader().loadTestsFromModule(tests);
    suite.addTests(thissuite)
#
unittest.TextTestRunner(verbosity=2).run(suite)
