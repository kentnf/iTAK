import argparse
import hashlib
import shutil
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import itak.runtime as itak_runtime


class ITAKRuntimeTests(unittest.TestCase):
    def test_build_analysis_request_converts_namespace(self):
        args = argparse.Namespace(
            seq_files=["a.fa", "b.fa"],
            frame="3R",
            process=4,
            specific="NAC",
            output="/tmp/out",
            mode="normal",
            classify=True,
        )

        request = itak_runtime.build_analysis_request(args)

        self.assertEqual(request.seq_files, ("a.fa", "b.fa"))
        self.assertEqual(request.frame, "3R")
        self.assertEqual(request.process, 4)
        self.assertEqual(request.specific, "NAC")
        self.assertEqual(request.output, "/tmp/out")
        self.assertEqual(request.mode, "normal")
        self.assertTrue(request.classify)

    def test_build_analysis_request_returns_existing_request(self):
        request = itak_runtime.AnalysisRequest(
            seq_files=("a.fa",),
            frame="6",
            process=1,
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
        self.assertIsNone(request.specific)
        self.assertIsNone(request.output)
        self.assertEqual(request.mode, "quick")
        self.assertFalse(request.classify)

    def test_find_database_source_dir_discovers_nested_database_folder(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            extract_root = Path(temp_dir)
            database_dir = extract_root / "custom-root" / "payload" / "database"
            database_dir.mkdir(parents=True)
            for file_name in itak_runtime.DATABASE_MARKER_FILES:
                (database_dir / file_name).write_text("ok\n")

            resolved = itak_runtime.find_database_source_dir(extract_root)

        self.assertEqual(resolved, database_dir)

    def test_resolve_database_location_skips_invalid_cwd_database_and_finds_conda_install(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            cwd_database_dir = temp_path / "workspace" / "database"
            cwd_database_dir.mkdir(parents=True)
            (cwd_database_dir / "db.txt").write_text("placeholder\n")

            conda_database_dir = temp_path / "pixi" / "share" / "itak" / "database"
            conda_database_dir.mkdir(parents=True)
            for file_name in itak_runtime.DATABASE_MARKER_FILES:
                (conda_database_dir / file_name).write_text("ok\n")

            with mock.patch("itak.runtime.os.getcwd", return_value=str(temp_path / "workspace")), \
                 mock.patch("itak.runtime.os.getenv", side_effect=lambda key: str(temp_path / "pixi") if key == "CONDA_PREFIX" else None):
                resolved = itak_runtime.resolve_database_location()

        self.assertEqual(resolved, itak_runtime.DatabaseLocation(path=str(conda_database_dir), source="conda"))

    def test_resolve_database_location_finds_user_install(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            user_database_dir = temp_path / ".local" / "share" / "itak" / "database"
            user_database_dir.mkdir(parents=True)
            for file_name in itak_runtime.DATABASE_MARKER_FILES:
                (user_database_dir / file_name).write_text("ok\n")

            with mock.patch("itak.runtime.get_user_database_dir", return_value=str(user_database_dir)), \
                 mock.patch("itak.runtime.os.getcwd", return_value=str(temp_path / "workspace")), \
                 mock.patch("itak.runtime.os.getenv", return_value=None):
                resolved = itak_runtime.resolve_database_location()

        self.assertEqual(resolved, itak_runtime.DatabaseLocation(path=str(user_database_dir), source="user"))

    def test_resolve_database_target_dir_matches_conda_default(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            conda_prefix = Path(temp_dir) / "pixi"

            def fake_getenv(key):
                if key == "CONDA_PREFIX":
                    return str(conda_prefix)
                return None

            with mock.patch("itak.runtime.os.getenv", side_effect=fake_getenv), \
                 mock.patch("itak.runtime.is_database_dir", return_value=False):
                resolved = itak_runtime.resolve_database_target_dir()

        self.assertEqual(resolved, str(conda_prefix / "share" / "itak" / "database"))

    def test_install_database_archive_accepts_custom_archive_root(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            source_root = temp_path / "custom-db-root" / "database"
            source_root.mkdir(parents=True)
            (source_root / "TF_Rule.txt").write_text("rule\n")
            (source_root / "GA_table.txt").write_text("ga\n")
            (source_root / "PK_class_desc.txt").write_text("pk\n")
            archive_path = temp_path / "itak-db-v1.tar.gz"
            with tarfile.open(archive_path, "w:gz") as archive_handle:
                archive_handle.add(temp_path / "custom-db-root", arcname="custom-db-root")

            sha256_path = temp_path / "itak-db-v1.tar.gz.sha256"
            sha256_path.write_text(f"{hashlib.sha256(archive_path.read_bytes()).hexdigest()}  {archive_path.name}\n")

            target_dir = temp_path / "target-db"

            def fake_download(url, destination_path):
                source_path = archive_path if url.endswith(".tar.gz") else sha256_path
                shutil.copy2(source_path, destination_path)

            with mock.patch("itak.runtime.download_file", side_effect=fake_download):
                itak_runtime.install_database_archive(
                    "https://example.com/itak-db-v1.tar.gz",
                    target_dir,
                    sha256_url="https://example.com/itak-db-v1.tar.gz.sha256",
                )

            self.assertEqual((target_dir / "TF_Rule.txt").read_text(), "rule\n")
            self.assertEqual((target_dir / "GA_table.txt").read_text(), "ga\n")
            self.assertEqual((target_dir / "PK_class_desc.txt").read_text(), "pk\n")


if __name__ == "__main__":
    unittest.main()
