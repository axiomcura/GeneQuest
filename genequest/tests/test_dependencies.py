# ----------------------------------------
# test_dependencies.py
# Author: Erik Serrano
# Email: erik.serrano@cuanschutz.edu
#
# checks if the correct dependencies and versions are installed properly
# ----------------------------------------
import unittest
import logging
from packaging import version


class TestDependencies(unittest.TestCase):
    """Testing all dependencies and versions installed"""

    # creating a stdout logger
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        format="%(asctime)s %(module)s %(levelname)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.INFO,
    )

    def test_numpy_import(self):
        """Checks if numpy and correct version is installed"""
        # import checking
        import_check = True
        try:
            import numpy as np
        except ImportError:
            import_check = False
            self.logger.error(
                "test_numpy_import import_check FAILED - unable to import numpy"
            )
            self.assertEqual(True, import_check)

        self.assertEqual(True, import_check)
        self.logger.info("test_numpy_import import_check: PASSED")

        # checking version
        req_version = version.parse("1.22.2")
        current_version = version.parse(np.__version__)
        version_check = current_version == req_version

        if version_check is False:
            self.logger.error(
                "test_numpy_import version_check: FAILED - incorrect version installed"
            )
            self.assertEqual(True, import_check)

        self.logger.info("test_numpy_import version_check: PASSED")
        self.assertEqual(True, version_check)

    def test_packaging_import(self):
        """Checks if packaging module and correct version is installed"""
        import_check = True
        try:
            import packaging
        except ImportError:
            import_check = False
            self.logger.error(
                "test_packaging_import: FAILED - unable to import packaging"
            )
            self.assertEqual(True, import_check)

        self.logger.info("test_packaging_import import_check: PASSED")
        self.assertEqual(True, import_check)

        # checking version
        req_version = version.parse("21.3")
        current_version = version.parse(packaging.__version__)
        version_check = current_version == req_version

        if version_check is False:
            self.logger.error(
                "test_numpy_import version_check: FAILED - incorrect version installed"
            )
            self.assertEqual(True, import_check)

        self.assertEqual(True, version_check)
        self.logger.info("test_numpy_import import_check: PASSED")

    def test_pandas_import(self):
        """Checks if pandas module and correct version is installed"""
        import_check = True
        try:
            import pandas as pd
        except ImportError:
            import_check = False
            self.logger.error("test_pandas_import: FAILED - unable to import pandas")
            self.assertEqual(True, import_check)

        # checking version
        req_version = version.parse("1.4.1")
        current_version = version.parse(pd.__version__)
        version_check = current_version == req_version

        if version_check is False:
            self.logger.error(
                "test_pandas_import version_check: FAILED - incorrect version installed"
            )
            self.assertEqual(True, import_check)

        self.assertEqual(True, version_check)
        self.logger.info("test_pandas_import import_check: PASSED")

    def test_memory_profiler_import(self):
        """Checks if pandas module and correct version is installed"""
        import_check = True
        try:
            import memory_profiler
        except ImportError:
            import_check = False
            self.logger.error(
                "test_memory_profiler_import: FAILED - unable to import memory_profiler"
            )
            self.assertEqual(True, import_check)

        # checking version
        req_version = version.parse("0.58.0")
        current_version = version.parse(memory_profiler.__version__)
        version_check = current_version == req_version

        if version_check is False:
            self.logger.error(
                "test_memory_profiler_import version_check: FAILED - incorrect version installed"
            )
            self.logger.info(
                "Please install memory_profiler by entering: pip install pip install memory-profiler=0.58.0"
            )
            self.assertEqual(True, import_check)

        self.assertEqual(True, version_check)
        self.logger.info("test_memory_profiler_import import_check: PASSED")
