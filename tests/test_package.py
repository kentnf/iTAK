import importlib
import unittest

import iTAK
import itak


class ITAKPackageTests(unittest.TestCase):
    def test_package_reexports_compatibility_api(self):
        self.assertEqual(itak.version, iTAK.version)
        self.assertEqual(itak.db_version, iTAK.db_version)
        self.assertIs(itak.translate_frames, iTAK.translate_frames)
        self.assertIs(itak.build_analysis_request, iTAK.build_analysis_request)

    def test_package_main_module_exposes_main(self):
        main_module = importlib.import_module("itak.__main__")

        self.assertTrue(callable(main_module.main))


if __name__ == "__main__":
    unittest.main()
