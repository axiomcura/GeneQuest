# ----------------------------------------
# run_test.py
# Author: Erik Serrano
# Email: erik.serrano@cuanschutz.edu
#
# main script for handing all testing components sequentially.
# steps:
# - Dependencies and version check
# - Testing function independently
# - Testing use cases
# ----------------------------------------
import unittest
import subprocess
import sys

# importing test modules
sys.path.append(".")
import test_dependencies
import test_functions

# creating test suite


def dependencies_test():
    """Tests Dependencies and versions"""
    # create testing env
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # adding test
    suite.addTest(loader.loadTestsFromModule(test_dependencies))

    # create a test runner and use test suite
    runner = unittest.TextTestRunner(verbosity=0)
    result = runner.run(suite)
    return result


def function_test():
    """Tests all functions developed in GeneQuest"""
    # create testing env
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # adding test
    suite.addTest(loader.loadTestsFromModule(test_functions))

    # create a test runner and use test suite
    runner = unittest.TextTestRunner(verbosity=0)
    result = runner.run(suite)
    return result


if __name__ == "__main__":

    # sequential test handler
    func_list = [function_test, dependencies_test]

    for test_func in func_list:
        results = test_func()
        if len(results.failures) > 0:

            print("here are the failed functions")
            for x in results.failures:
                failed_func = str(x[0]).split()[0]
                print(failed_func)
            exit()
