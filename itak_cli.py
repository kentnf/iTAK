import tempfile
from pathlib import Path

from itak.cli import DB_URL_LATEST, DATABASE_FILES, db_version, debug, send_mail, version
from itak.pipeline import itak_identify
from itak.runtime import build_arg_parser, download_and_extract, prepare_runtime_environment


def run_cli(args):
    bin_path, dbs_path = prepare_runtime_environment(DATABASE_FILES, DB_URL_LATEST)

    if args.update:
        with tempfile.TemporaryDirectory() as temp_dir:
            download_and_extract(DB_URL_LATEST, Path(temp_dir), dbs_path)
    else:
        itak_identify(args, bin_path, dbs_path, debug=debug)


def main():
    args = build_arg_parser().parse_args()
    run_cli(args)
