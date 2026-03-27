import argparse
import unittest
from unittest import mock

import itak_cli


class ITAKCliTests(unittest.TestCase):
    def test_main_parses_args_before_runtime_setup(self):
        parser = mock.Mock()
        parser.parse_args.side_effect = SystemExit(0)

        with mock.patch.object(itak_cli, "build_arg_parser", return_value=parser), \
             mock.patch.object(itak_cli, "prepare_runtime_environment") as prepare_runtime_environment:
            with self.assertRaises(SystemExit):
                itak_cli.main()

        prepare_runtime_environment.assert_not_called()

    def test_run_cli_dispatches_identify_path(self):
        args = argparse.Namespace(update=False)

        with mock.patch.object(itak_cli, "prepare_runtime_environment", return_value=("bin", "db")) as prepare_runtime_environment, \
             mock.patch.object(itak_cli, "itak_identify") as itak_identify:
            itak_cli.run_cli(args)

        prepare_runtime_environment.assert_called_once_with(itak_cli.DATABASE_FILES, itak_cli.DB_URL_LATEST)
        itak_identify.assert_called_once_with(args, "bin", "db", debug=itak_cli.debug)

    def test_run_cli_dispatches_update_path(self):
        args = argparse.Namespace(update=True)

        with mock.patch.object(itak_cli, "prepare_runtime_environment", return_value=("bin", "db")) as prepare_runtime_environment, \
             mock.patch.object(itak_cli, "download_and_extract") as download_and_extract:
            itak_cli.run_cli(args)

        prepare_runtime_environment.assert_called_once_with(itak_cli.DATABASE_FILES, itak_cli.DB_URL_LATEST)
        download_and_extract.assert_called_once()
        self.assertEqual(download_and_extract.call_args.args[0], itak_cli.DB_URL_LATEST)
        self.assertEqual(download_and_extract.call_args.args[2], "db")


if __name__ == "__main__":
    unittest.main()
