import os
import smtplib
import tempfile
from email.mime.text import MIMEText
from pathlib import Path

from itak.pipeline import itak_identify
from itak.runtime import build_arg_parser, download_and_extract, prepare_runtime_environment


version = '2.0.2'
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

DB_URL_LATEST = "https://github.com/kentnf/iTAK/archive/refs/tags/db-v1.tar.gz"


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
    bin_path, dbs_path = prepare_runtime_environment(DATABASE_FILES, DB_URL_LATEST)

    if args.update:
        with tempfile.TemporaryDirectory() as temp_dir:
            download_and_extract(DB_URL_LATEST, Path(temp_dir), dbs_path)
    else:
        itak_identify(args, bin_path, dbs_path, debug=debug)


def main():
    args = build_arg_parser().parse_args()
    run_cli(args)
