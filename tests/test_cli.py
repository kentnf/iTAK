import argparse
import io
import unittest
from contextlib import redirect_stdout
from unittest import mock

import itak.cli as itak_cli


class ITAKCliTests(unittest.TestCase):
    def test_normalize_argv_strips_leading_separator(self):
        self.assertEqual(itak_cli.normalize_argv(["--", "db", "path"]), ["db", "path"])
        self.assertEqual(itak_cli.normalize_argv(["test_seq"]), ["test_seq"])

    def test_main_parses_args_before_runtime_setup(self):
        parser = mock.Mock()
        parser.parse_args.side_effect = SystemExit(0)

        with mock.patch("itak.runtime.build_arg_parser", return_value=parser), \
             mock.patch.object(itak_cli, "prepare_runtime_environment") as prepare_runtime_environment:
            with self.assertRaises(SystemExit):
                itak_cli.main()

        prepare_runtime_environment.assert_not_called()

    def test_main_accepts_leading_separator_for_db_commands(self):
        args = argparse.Namespace(db_command="path", target=True)

        with mock.patch.object(itak_cli, "build_db_parser") as build_db_parser, \
             mock.patch.object(itak_cli, "run_db_cli") as run_db_cli, \
             mock.patch("itak.cli.sys.argv", ["itak", "--", "db", "path", "--target"]):
            build_db_parser.return_value.parse_args.return_value = args
            itak_cli.main()

        build_db_parser.return_value.parse_args.assert_called_once_with(["path", "--target"])
        run_db_cli.assert_called_once_with(args)

    def test_run_cli_dispatches_identify_path(self):
        args = argparse.Namespace()

        with mock.patch.object(itak_cli, "prepare_runtime_environment", return_value=("bin", "db")) as prepare_runtime_environment, \
             mock.patch.object(itak_cli, "itak_identify") as itak_identify:
            itak_cli.run_cli(args)

        prepare_runtime_environment.assert_called_once_with(itak_cli.DATABASE_FILES)
        itak_identify.assert_called_once_with(args, "bin", "db", debug=itak_cli.debug)

    def test_run_db_cli_path_prints_resolved_location(self):
        args = argparse.Namespace(db_command="path", target=False)

        with mock.patch("itak.cli.resolve_database_location", return_value=argparse.Namespace(path="/tmp/itak-db", source="env")):
            output = io.StringIO()
            with redirect_stdout(output):
                itak_cli.run_db_cli(args)

        self.assertEqual(output.getvalue().strip(), "/tmp/itak-db")

    def test_run_db_cli_verify_prefers_resolved_location(self):
        args = argparse.Namespace(db_command="verify", db_path=None)

        with mock.patch("itak.cli.resolve_database_location", return_value=argparse.Namespace(path="/resolved/db", source="env")), \
             mock.patch("itak.cli.resolve_database_target_dir", return_value="/target/db") as resolve_database_target_dir, \
             mock.patch("itak.cli.verify_database_files", return_value={}) as verify_database_files:
            output = io.StringIO()
            with redirect_stdout(output):
                itak_cli.run_db_cli(args)

        resolve_database_target_dir.assert_not_called()
        verify_database_files.assert_called_once_with("/resolved/db", itak_cli.DATABASE_FILES)
        self.assertIn("Database verification passed: /resolved/db", output.getvalue())


if __name__ == "__main__":
    unittest.main()
