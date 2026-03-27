import argparse
import unittest

import itak_runtime


class ITAKRuntimeTests(unittest.TestCase):
    def test_build_analysis_request_converts_namespace(self):
        args = argparse.Namespace(
            seq_files=["a.fa", "b.fa"],
            frame="3R",
            process=4,
            update=True,
            specific="NAC",
            output="/tmp/out",
            mode="normal",
            classify=True,
        )

        request = itak_runtime.build_analysis_request(args)

        self.assertEqual(request.seq_files, ("a.fa", "b.fa"))
        self.assertEqual(request.frame, "3R")
        self.assertEqual(request.process, 4)
        self.assertTrue(request.update)
        self.assertEqual(request.specific, "NAC")
        self.assertEqual(request.output, "/tmp/out")
        self.assertEqual(request.mode, "normal")
        self.assertTrue(request.classify)

    def test_build_analysis_request_returns_existing_request(self):
        request = itak_runtime.AnalysisRequest(
            seq_files=("a.fa",),
            frame="6",
            process=1,
            update=False,
            specific=None,
            output=None,
            mode="quick",
            classify=False,
        )

        self.assertIs(itak_runtime.build_analysis_request(request), request)

    def test_build_analysis_request_fills_defaults_for_missing_fields(self):
        args = argparse.Namespace(seq_files=["a.fa"])

        request = itak_runtime.build_analysis_request(args)

        self.assertEqual(request.seq_files, ("a.fa",))
        self.assertEqual(request.frame, "6")
        self.assertEqual(request.process, 1)
        self.assertFalse(request.update)
        self.assertIsNone(request.specific)
        self.assertIsNone(request.output)
        self.assertEqual(request.mode, "quick")
        self.assertFalse(request.classify)


if __name__ == "__main__":
    unittest.main()
