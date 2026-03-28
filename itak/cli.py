import argparse
import os
import smtplib
import sys
from email.mime.text import MIMEText

from itak import __version__
from itak.pipeline import itak_identify
from itak.runtime import install_database_archive, prepare_runtime_environment, resolve_database_location, resolve_database_target_dir, verify_database_files


version = __version__
db_version = '2.1'
debug = False

DATABASE_FILES = {
    'rule': ["TF_Rule.txt", "GA_table.txt", "PK_class_desc.txt"],
    'hmm': [
        "Tfam_domain.hmm",
        "TF_selfbuild.hmm",
        "PlantsPHMM3_89.hmm",
        "Plant_Pkinase_fam.hmm",
        "Pkinase_sub_WNK1.hmm",
        "Pkinase_sub_MAK.hmm",
    ],
    'other': ["TF_selfbuild.sto"],
}

DB_URL_LATEST = "https://github.com/kentnf/iTAK/releases/download/db-v1/itak-db-v1.tar.gz"
DB_SHA256_URL_LATEST = "https://github.com/kentnf/iTAK/releases/download/db-v1/itak-db-v1.tar.gz.sha256"


def send_mail(address, input_file):
    job_id = os.path.basename(input_file)
    download_link = f"http://bioinfo.bti.cornell.edu/cgi-bin/itak/online_itak.cgi?rid={job_id}"
    message = f"Hi,\n the analysis for {job_id} is finished. Please view and download your result through link {download_link}\nThank you for using iTAK.\n"

    msg = MIMEText(message)
    msg['To'] = address
    msg['From'] = 'bioinfo@cornell.edu'
    msg['Subject'] = f"[iTAK] analysis for {job_id} is finished"

    try:
        server = smtplib.SMTP('localhost')
        server.sendmail('bioinfo@cornell.edu', [address], msg.as_string())
        server.quit()
        print("Mail sent successfully.")
    except Exception as e:
        print(f"Failed to send mail: {e}")


def run_cli(args):
    bin_path, dbs_path = prepare_runtime_environment(DATABASE_FILES)
    itak_identify(args, bin_path, dbs_path, debug=debug)


def build_db_parser():
    parser = argparse.ArgumentParser(description="iTAK database management")
    db_subparsers = parser.add_subparsers(dest="db_command", required=True)

    db_path_parser = db_subparsers.add_parser("path", help="Print the resolved database path")
    db_path_parser.add_argument("--target", action="store_true", help="Print the default install target path")

    db_verify_parser = db_subparsers.add_parser("verify", help="Verify required database files")
    db_verify_parser.add_argument("--path", dest="db_path", help="Database path to verify")

    db_download_parser = db_subparsers.add_parser("download", help="Download and install database files")
    db_download_parser.add_argument("--path", dest="db_path", help="Install database into this directory")
    db_download_parser.add_argument("--url", default=DB_URL_LATEST, help="Database archive URL")
    db_download_parser.add_argument("--sha256-url", default=DB_SHA256_URL_LATEST, help="SHA256 file URL")

    return parser


def run_db_cli(args):
    if args.db_command == "path":
        if args.target:
            print(resolve_database_target_dir())
            return
        location = resolve_database_location()
        if location is None:
            print("Database not found")
            return
        print(location.path)
        return

    if args.db_command == "verify":
        if args.db_path:
            db_path = args.db_path
        else:
            location = resolve_database_location()
            db_path = location.path if location is not None else resolve_database_target_dir()
        missing = verify_database_files(db_path, DATABASE_FILES)
        if missing:
            print(f"Database verification failed for: {db_path}")
            for file_type in sorted(missing.keys()):
                print(f"{file_type}: {', '.join(missing[file_type])}")
            raise SystemExit(1)
        print(f"Database verification passed: {db_path}")
        return

    if args.db_command == "download":
        db_path = resolve_database_target_dir(args.db_path)
        install_database_archive(args.url, db_path, sha256_url=args.sha256_url)
        print(f"Database installed to: {db_path}")
        return

    raise SystemExit(f"Unknown db command: {args.db_command}")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "db":
        args = build_db_parser().parse_args(sys.argv[2:])
        run_db_cli(args)
        return
    from itak.runtime import build_arg_parser

    args = build_arg_parser().parse_args()
    run_cli(args)
