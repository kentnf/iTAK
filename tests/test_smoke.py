import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class ITAKSmokeTest(unittest.TestCase):
    def test_cli_matches_expected_fixture_output(self):
        repo_root = Path(__file__).resolve().parents[1]
        script_path = repo_root / "iTAK.py"
        input_path = repo_root / "test_seq"
        expected_dir = repo_root / "test_seq_output"

        required_files = [
            "pk_sequence.fasta",
            "shiu_alignment.txt",
            "shiu_classification.txt",
            "tf_alignment.txt",
            "tf_classification.txt",
            "tf_sequence.fasta",
        ]

        if shutil.which("hmmscan") is None:
            self.skipTest("hmmscan is required for the smoke test")

        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "output"

            subprocess.run(
                [sys.executable, str(script_path), str(input_path), "-o", str(output_dir)],
                cwd=repo_root,
                check=True,
            )

            for filename in required_files:
                with self.subTest(filename=filename):
                    expected_path = expected_dir / filename
                    actual_path = output_dir / filename

                    self.assertTrue(actual_path.exists(), f"missing output file: {filename}")
                    self.assertEqual(
                        actual_path.read_text(),
                        expected_path.read_text(),
                        f"output mismatch for {filename}",
                    )


if __name__ == "__main__":
    unittest.main()
