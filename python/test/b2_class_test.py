import mpfr_test
import function_tree_test
import system_test
import unittest














if __name__ == '__main__':
    mods = (mpfr_test,function_tree_test,system_test)
    suite = unittest.TestSuite();
    for tests in mods:
        thissuite = unittest.TestLoader().loadTestsFromModule(tests);
        suite.addTests(thissuite)

    unittest.TextTestRunner(verbosity=2).run(suite)