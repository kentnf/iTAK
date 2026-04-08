import importlib
import unittest

import itak
from itak.runtime import build_analysis_request
from itak.sequences import translate_frames


class ITAKPackageTests(unittest.TestCase):
    def test_package_exposes_version_metadata(self):
        self.assertEqual(itak.__version__, "2.0.6")
        self.assertTrue(callable(translate_frames))
        self.assertTrue(callable(build_analysis_request))

    def test_package_main_module_exposes_main(self):
        main_module = importlib.import_module("itak.__main__")

        self.assertTrue(callable(main_module.main))


if __name__ == "__main__":
    unittest.main()
