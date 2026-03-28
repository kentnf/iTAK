import argparse
import hashlib
import os
import platform
import requests
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

from itak.rules import load_rule_records

DEFAULT_DB_ENV_VAR = "ITAK_DB_DIR"


@dataclass(frozen=True)
class RuntimeConfig:
    mode: str
    frame: str
    cpu: str
    bin_dir: str
    dbs_dir: str
    hmmscan_bin: str
    tfam_db: str
    sfam_db: str
    plantsp_db: str
    shiu_db: str
    pfam_db: str
    correct_ga: str
    tf_rule: dict
    ga_cutoff: dict
    pkid_des: dict


@dataclass(frozen=True)
class SamplePaths:
    input_file: str
    temp_dir: str
    output_dir: str
    input_protein_f: str
    tmp_pfam_hmmscan: str


@dataclass(frozen=True)
class AnalysisRequest:
    seq_files: tuple[str, ...]
    frame: str
    process: int
    specific: str | None
    output: str | None
    mode: str
    classify: bool


@dataclass(frozen=True)
class DatabaseLocation:
    path: str
    source: str


def check_os():
    try:
        os_name = platform.system()
        if os_name not in ["Linux", "Darwin"]:
            raise Exception(f"Unsupported operating system: {os_name}")
        print(f"Running on supported OS: {os_name}")
    except Exception as e:
        print(f"Error: {e}")
        print(f"iTAK does not support running on {os_name}")
        sys.exit(1)


def check_hmmer():
    try:
        subprocess.run(["hmmscan", "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        print("HMMER is installed.")
        return os.path.dirname(shutil.which("hmmscan"))
    except subprocess.CalledProcessError:
        print("HMMER is installed, but there was a problem running it.")
        sys.exit(1)
    except FileNotFoundError:
        print("HMMER is not installed.")
        sys.exit(1)


def check_db_file_exit(db_path, database_files, ftype, bin_path):
    missing_files = []
    for file_name in database_files[ftype]:
        full_path = os.path.join(db_path, file_name)
        if os.path.exists(full_path):
            if ftype == 'hmm':
                if not all(os.path.isfile(f"{full_path}.{ext}") for ext in ["h3f", "h3i", "h3m", "h3p"]):
                    print(f"[WARN]no database file {full_path}")
                    subprocess.run([os.path.join(bin_path, "hmmpress"), "-f", full_path], check=True)
                else:
                    print(f"Database file {full_path} is ready.")
        else:
            missing_files.append(file_name)
            print(f"File does not exist: {full_path}")

    return missing_files


def download_file(url, destination_path):
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        with open(destination_path, 'wb') as output_file:
            for chunk in response.iter_content(chunk_size=8192):
                output_file.write(chunk)


def compute_sha256(file_path):
    sha256 = hashlib.sha256()
    with open(file_path, 'rb') as file_handle:
        for chunk in iter(lambda: file_handle.read(1024 * 1024), b''):
            sha256.update(chunk)
    return sha256.hexdigest()


def read_expected_sha256(sha256_file):
    with open(sha256_file) as file_handle:
        first_line = file_handle.readline().strip()
    return first_line.split()[0]


def check_specific_families(dbs_path):
    spec_path = dbs_path + '/specific'
    spec_families = []
    for hmm_file in [file_name for file_name in os.listdir(spec_path) if file_name.endswith('.hmm')]:
        txt_path = os.path.join(spec_path, os.path.splitext(hmm_file)[0] + '.txt')
        if os.path.exists(txt_path):
            spec_families.append(os.path.splitext(hmm_file)[0])
    return spec_families


def run_cmd(cmd, debug=False):
    if not cmd:
        print("[ERR]no command")
        exit()

    if debug:
        print(cmd)
        return True

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        print(f"[ERR]cmd: {cmd}")
        raise


def get_db_path():
    location = resolve_database_location()
    if location is not None:
        return location.path

    raise FileNotFoundError(
        "Can not find the database path. Run `itak db download` to install the database files "
        f"or set {DEFAULT_DB_ENV_VAR}."
    )


def resolve_database_location():
    env_db_path = os.getenv(DEFAULT_DB_ENV_VAR)
    if env_db_path and os.path.exists(env_db_path):
        return DatabaseLocation(path=env_db_path, source=DEFAULT_DB_ENV_VAR)

    current_dir_dbpath = os.path.join(os.getcwd(), 'database')
    if os.path.exists(current_dir_dbpath):
        return DatabaseLocation(path=current_dir_dbpath, source="cwd")

    conda_prefix = os.getenv('CONDA_PREFIX')
    if conda_prefix:
        db_path = os.path.join(conda_prefix, 'share', 'itak', 'database')
        if os.path.exists(db_path):
            return DatabaseLocation(path=db_path, source="conda")

    return None


def resolve_database_target_dir(custom_path=None):
    if custom_path is not None:
        return custom_path

    env_db_path = os.getenv(DEFAULT_DB_ENV_VAR)
    if env_db_path:
        return env_db_path

    conda_prefix = os.getenv('CONDA_PREFIX')
    if conda_prefix:
        return os.path.join(conda_prefix, 'share', 'itak', 'database')

    return os.path.join(Path.home(), '.local', 'share', 'itak', 'database')


def find_database_source_dir(extract_root):
    extract_root = Path(extract_root)
    candidates = [extract_root]
    candidates.extend(path for path in extract_root.rglob("database") if path.is_dir())

    for candidate in candidates:
        if all((candidate / file_name).exists() for file_name in ("TF_Rule.txt", "GA_table.txt", "PK_class_desc.txt")):
            return candidate

    raise FileNotFoundError("Can not find a database directory in the extracted archive")


def file_has_content(path):
    return os.path.exists(path) and os.path.getsize(path) > 0


def resolve_output_dir(input_file, custom_output, multiple_inputs):
    if custom_output is None:
        return f"{input_file}_output"
    if multiple_inputs:
        return os.path.join(custom_output, f"{Path(input_file).name}_output")
    return custom_output


def build_hmmscan_cmd(hmmscan_bin, cpu, output_path, db_path, input_path):
    return [
        hmmscan_bin,
        "--acc",
        "--notextw",
        "--cpu",
        str(cpu),
        "-o",
        output_path,
        db_path,
        input_path,
    ]


def resolve_cpu(process_count):
    if process_count > 0:
        return str(process_count)
    return '20'


def resolve_frame_mode(frame_arg):
    if frame_arg in ('3F', '3R'):
        return frame_arg
    return '6'


def validate_input_files(seq_files, custom_output):
    multiple_inputs = len(seq_files) > 1
    for input_file in seq_files:
        output_dir = resolve_output_dir(input_file, custom_output, multiple_inputs)
        temp_dir = f"{input_file}_temp"
        if os.path.exists(output_dir):
            print(f"[WARN]output folder exist: {output_dir}")
        if os.path.exists(temp_dir):
            print(f"[WARN]temp folder exist: {temp_dir}")
        if not os.path.isfile(input_file) or os.path.getsize(input_file) == 0:
            print("[ERR]input file not exist")
            sys.exit(1)


def pk_to_hash(file_path):
    hash_ = {}
    try:
        with open(file_path, 'r') as pfh:
            for line in pfh:
                pm = line.strip().split('\t', 2)
                hash_[pm[0]] = pm[1]
    except IOError as e:
        raise IOError(f"Can not open protein kinase description file: {file_path} {e}")
    return hash_


def build_analysis_request(args):
    if isinstance(args, AnalysisRequest):
        return args

    return AnalysisRequest(
        seq_files=tuple(getattr(args, "seq_files", ()) or ()),
        frame=getattr(args, "frame", "6"),
        process=getattr(args, "process", 1),
        specific=getattr(args, "specific", None),
        output=getattr(args, "output", None),
        mode=getattr(args, "mode", "quick"),
        classify=bool(getattr(args, "classify", False)),
    )


def build_runtime_config(args, bin_dir, dbs_dir, load_ga_cutoff_fn):
    request = build_analysis_request(args)
    mode = request.mode
    tfam_db = os.path.join(dbs_dir, "Tfam_domain.hmm")
    sfam_db = os.path.join(dbs_dir, "TF_selfbuild.hmm")
    plantsp_db = os.path.join(dbs_dir, "PlantsPHMM3_89.hmm")
    shiu_db = os.path.join(dbs_dir, "Plant_Pkinase_fam.hmm")
    pfam_db = os.path.join(dbs_dir, "Pfam-A.hmm")
    correct_ga = os.path.join(dbs_dir, "GA_table.txt")
    pk_desc = os.path.join(dbs_dir, "PK_class_desc.txt")

    tf_rule = load_rule_records(os.path.join(dbs_dir, "TF_Rule.txt"))
    if mode == 'normal':
        ga_cutoff = load_ga_cutoff_fn(pfam_db, correct_ga, sfam_db)
    elif mode == 'quick':
        ga_cutoff = load_ga_cutoff_fn(tfam_db, correct_ga, sfam_db)
    else:
        print(f"Mode error: {mode}")
        sys.exit(1)

    return RuntimeConfig(
        mode=mode,
        frame=resolve_frame_mode(request.frame),
        cpu=resolve_cpu(request.process),
        bin_dir=bin_dir,
        dbs_dir=dbs_dir,
        hmmscan_bin=os.path.join(bin_dir, "hmmscan"),
        tfam_db=tfam_db,
        sfam_db=sfam_db,
        plantsp_db=plantsp_db,
        shiu_db=shiu_db,
        pfam_db=pfam_db,
        correct_ga=correct_ga,
        tf_rule=tf_rule,
        ga_cutoff=ga_cutoff,
        pkid_des=pk_to_hash(pk_desc),
    )


def build_sample_paths(input_file, custom_output, multiple_inputs):
    temp_dir = f"{input_file}_temp"
    output_dir = resolve_output_dir(input_file, custom_output, multiple_inputs)
    return SamplePaths(
        input_file=input_file,
        temp_dir=temp_dir,
        output_dir=output_dir,
        input_protein_f=os.path.join(temp_dir, "protein_seq.fa"),
        tmp_pfam_hmmscan=os.path.join(temp_dir, "protein_seq.pfam.hmmscan.txt"),
    )


def ensure_sample_dirs(paths):
    os.makedirs(paths.temp_dir, exist_ok=True)
    os.makedirs(paths.output_dir, exist_ok=True)


def prepare_runtime_environment(database_files):
    check_os()
    bin_path = check_hmmer()
    try:
        dbs_path = get_db_path()
    except FileNotFoundError as exc:
        print(exc)
        sys.exit(1)

    missing_hmm = check_db_file_exit(dbs_path, database_files, "hmm", bin_path)
    if len(missing_hmm) > 0:
        print(f"Can not find database files: {missing_hmm}")
        print("Run `itak db download` to install the database files.")
        sys.exit(1)

    missing_rule = check_db_file_exit(dbs_path, database_files, "rule", bin_path)
    if len(missing_rule) > 0:
        print(f"Can not find rule file: {missing_rule}")
        sys.exit(1)

    return bin_path, dbs_path


def verify_database_files(dbs_path, database_files):
    missing = {}
    for file_type in ("rule", "hmm"):
        missing_files = []
        for file_name in database_files[file_type]:
            full_path = os.path.join(dbs_path, file_name)
            if not os.path.exists(full_path):
                missing_files.append(file_name)
        if missing_files:
            missing[file_type] = missing_files
    return missing


def install_database_archive(url, dbs_path, sha256_url=None):
    target_dir = Path(dbs_path)
    target_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        archive_path = temp_dir / url.split('/')[-1]
        download_file(url, archive_path)

        if sha256_url is not None:
            sha256_path = temp_dir / Path(sha256_url).name
            download_file(sha256_url, sha256_path)
            expected_sha256 = read_expected_sha256(sha256_path)
            actual_sha256 = compute_sha256(archive_path)
            if actual_sha256 != expected_sha256:
                raise ValueError(f"SHA256 mismatch for {archive_path.name}")

        subprocess.run(['tar', '-xzf', archive_path, '-C', temp_dir], check=True)
        source_dir = find_database_source_dir(temp_dir)
        for source_path in source_dir.iterdir():
            target_path = target_dir / source_path.name
            if source_path.is_dir():
                shutil.copytree(source_path, target_path, dirs_exist_ok=True)
            else:
                shutil.copy2(source_path, target_path)


def build_arg_parser():
    parser = argparse.ArgumentParser(description='iTAK -- Plant Transcription factor & Protein Kinase Identifier and Classifier')
    parser.add_argument('seq_files', nargs='+', help='Input sequence file(s)')
    parser.add_argument('-f', '--frame', help='Translate frame. (3F, 3R, 6; default = 6)', default='6')
    parser.add_argument('-p', '--process', type=int, help='Number of CPUs used for hmmscan. (default = 1)', default=1)
    parser.add_argument('-s', '--specific', help='User-defined specific family name for identification and classification')
    parser.add_argument('-o', '--output', help='Name of the output directory. (default = \'input file name\' + \'_output\')')
    parser.add_argument('-m', '--mode', help='Mode, quick or normal, please do not change it. (default = quick)', default='quick')
    parser.add_argument('-c', '--classify', action='store_true', help='OBSOLETE: Enable protein kinase PPC classification. (default = disable)')
    return parser
