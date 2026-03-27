import sys
import tempfile

from itak.cli import (
    DB_SHA256_URL_LATEST,
    DB_URL_LATEST,
    DATABASE_FILES,
    build_db_parser,
    db_version,
    debug,
    run_db_cli,
    send_mail,
    version,
)
from itak.pipeline import itak_identify
from itak.runtime import build_arg_parser, install_database_archive, prepare_runtime_environment


def run_cli(args):
    bin_path, dbs_path = prepare_runtime_environment(DATABASE_FILES, DB_URL_LATEST)

    if args.update:
        with tempfile.TemporaryDirectory() as temp_dir:
            install_database_archive(DB_URL_LATEST, dbs_path, sha256_url=None)
    else:
        itak_identify(args, bin_path, dbs_path, debug=debug)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "db":
        args = build_db_parser().parse_args(sys.argv[2:])
        run_db_cli(args)
        return
    args = build_arg_parser().parse_args()
    run_cli(args)
