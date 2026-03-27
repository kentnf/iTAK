import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO


class ITAKSmokeTest(unittest.TestCase):
    def _read_first_column_set(self, path):
        values = set()
        with open(path) as file_handle:
            for line in file_handle:
                line = line.strip()
                if not line:
                    continue
                values.add(line.split("\t", 1)[0])
        return values

    def test_cli_matches_expected_fixture_output(self):
        repo_root = Path(__file__).resolve().parents[1]
        script_path = repo_root / "iTAK.py"
        fixture_input_path = repo_root / "test_seq"
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
            temp_dir_path = Path(temp_dir)
            input_path = temp_dir_path / fixture_input_path.name
            output_dir = Path(temp_dir) / "output"
            temp_work_dir = Path(f"{input_path}_temp")
            shutil.copy2(fixture_input_path, input_path)
            shutil.rmtree(temp_work_dir, ignore_errors=True)

            try:
                subprocess.run(
                    [sys.executable, str(script_path), str(input_path), "-o", str(output_dir)],
                    cwd=repo_root,
                    check=True,
                )
            finally:
                shutil.rmtree(temp_work_dir, ignore_errors=True)

            for filename in required_files:
                with self.subTest(filename=filename):
                    expected_path = expected_dir / filename
                    actual_path = output_dir / filename

                    self.assertTrue(actual_path.exists(), f"missing output file: {filename}")

                    if filename.startswith("tf_"):
                        self.assertEqual(
                            actual_path.read_text(),
                            expected_path.read_text(),
                            f"output mismatch for {filename}",
                        )
                        continue

                    if filename == "pk_sequence.fasta":
                        with open(expected_path) as expected_handle:
                            expected_ids = [record.id for record in SeqIO.parse(expected_handle, "fasta")]
                        with open(actual_path) as actual_handle:
                            actual_ids = [record.id for record in SeqIO.parse(actual_handle, "fasta")]
                        self.assertEqual(actual_ids, expected_ids, f"record ids mismatch for {filename}")
                        continue

                    expected_ids = self._read_first_column_set(expected_path)
                    actual_ids = self._read_first_column_set(actual_path)
                    self.assertEqual(actual_ids, expected_ids, f"record ids mismatch for {filename}")


if __name__ == "__main__":
    unittest.main()
