#!/usr/bin/env python

import sys, os, re, platform, argparse, tempfile, requests, subprocess, shutil, smtplib
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
# from subprocess import run, call, CalledProcessError
# from shutil import copytree, rmtree, which
from email.mime.text import MIMEText
from email.header import Header

import itakm

# Define the iTAK version and debug mode
version = '2.0.2'
db_version = '2.1' # db-v1 iTAK2 database version 1
debug = False


def check_package_installed(package_name):
    result = subprocess.run(['conda', 'list'], capture_output=True, text=True)
    installed_packages = [line.split()[0] for line in result.stdout.splitlines()[2:]]
    if package_name in installed_packages:
        return True
    else:
        return False


def check_os():
    try:
        os_name = platform.system()
        if os_name not in ["Linux", "Darwin"]:
            raise Exception(f"Unsupported operating system: {os_name}")
        else:
            print(f"Running on supported OS: {os_name}")
    except Exception as e:
        print(f"Error: {e}")
        print(f"iTAK does not support running on {os_name}")
        sys.exit(1)


def check_hmmer():
    try:
        # run hmmscan command 
        subprocess.run(["hmmscan", "-h"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        print("HMMER is installed.")
        hmmscan_path = os.path.dirname(shutil.which("hmmscan"))
        return hmmscan_path
    except subprocess.CalledProcessError:
        # if failed, return the state
        print("HMMER is installed, but there was a problem running it.")
        sys.exit(1)
    except FileNotFoundError:
        # if can not find hmmer
        print("HMMER is not installed.")
        sys.exit(1)

# the ftype could be: rule, hmm, other
def check_db_file_exit(db_path, database_files, ftype, bin_path):
    missing_files = []
    for file in database_files[ftype]:
        full_path = os.path.join(db_path, file)
        if os.path.exists(full_path):
            if ftype == 'hmm':
                # building database using hmmpress 
                if not all(os.path.isfile(f"{full_path}.{ext}") for ext in ["h3f", "h3i", "h3m", "h3p"]):
                    print(f"[WARN]no database file {full_path}")
                    subprocess.call([os.path.join(bin_path, "hmmpress"), "-f", full_path])
                else:
                    print(f"Database file {full_path} is ready.")
        else:
            missing_files.append(file)
            print(f"File does not exist: {full_path}")

    return missing_files


def download_and_extract(url, temp_dir, dbs_path):
    # download file
    local_filename = url.split('/')[-1]
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(temp_dir / local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    
    # uncompress the tar.gz and move files to database  
    subprocess.run(['tar', '-xzf', temp_dir / local_filename, '-C', temp_dir], check=True)
    subprocess.run(f'cp -rf {temp_dir}/iTAK-db-v1/database/* {dbs_path}', shell=True, check=True)

def check_and_prepare_database(dbs_path, url, database_files, bin_path):
    #and check/build db files
    missing_files = check_db_file_exit(dbs_path, database_files, 'hmm', bin_path)

    # dowload database files and check/build db files 
    if len(missing_files) > 0:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)
        
            # download and extract database files
            download_and_extract(url, temp_dir, dbs_path)
        
            # list of database files (for debug)
            #result = subprocess.run(['find', temp_dir, '-type', 'f'], capture_output=True, text=True)
            #if result.returncode == 0:
            #    print(result.stdout)
            #else:
            #    print("Error:", result.stderr)
            check_db_file_exit(dbs_path, database_files, 'hmm', bin_path)

def check_specific_families(dbs_path):
    spec_path = dbs_path + '/specific'
    hmm_files = [f for f in os.listdir(spec_path) if f.endswith('.hmm')]
    spec_families = []

    for hmm_file in hmm_files:
        txt_file = os.path.splitext(hmm_file)[0] + '.txt'
        txt_path = os.path.join(spec_path, txt_file)

        if os.path.exists(txt_path):
            spec_families.append(os.path.splitext(hmm_file)[0])

    return spec_families

# send mail to user when analysis is done, deprecated in standalone version 
# Example: send_mail("example@email.com", "/path/to/input_file")

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


# example: run_cmd("echo 'Hello World'", debug=True)
def run_cmd(cmd, debug=False):
    if not cmd:
        print("[ERR]no command")
        exit()

    if debug:
        print(cmd)
        return True

    try:
        subprocess.call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERR]cmd: {cmd}")
        raise

def get_db_path():
    # check conda prefix and db path in conda env
    conda_prefix = os.getenv('CONDA_PREFIX')
    # itak_installed = check_package_installed('itak')
    if conda_prefix:
        db_path = os.path.join(conda_prefix, 'share', 'itak', 'database')
        if os.path.exists(db_path):
            return db_path
        else:
            os.makedirs(db_path)
            return db_path

    # check current dir db path
    current_dir_dbpath = os.path.join(os.getcwd(), 'database')
    if os.path.exists(current_dir_dbpath):
        return current_dir_dbpath
    # can not find db path, reinstall the program
    print(f"Can not find the database path, please run iTAK.py to install database")
    #sys.exit(1)


def get_bin_path():
    # check conda prefix and bin path in conda env
    conda_prefix = os.getenv('CONDA_PREFIX')
    # itak_installed = check_package_installed('itak')
    # if conda_prefix and itak_installed:
    if conda_prefix:
        bin_path = os.path.join(conda_prefix, 'bin')
        if os.path.exists(bin_path):
            return bin_path
    # check current dir bin path
    current_dir_binpath = os.path.join(os.getcwd(), 'bin')
    if os.path.exists(current_dir_binpath):
        return current_dir_binpath

    # can not find bin path
    print(f"Can not find the bin path for hmmer, please re-install iTAK or install hmmer")
    sys.exit(1)


def is_nucleotide(sequence):
    nucleotide_set = set("ACGTUacgtu")
    return set(sequence).issubset(nucleotide_set)

# fasta to hash, key: id, value->alphabet: , value->seq: seq;
def seq_to_hash(input_file):
    """
    put seq info to hash
    """
    seq_info = {}
    in_ = SeqIO.parse(input_file, 'fasta')
    for inseq in in_:
        alphabet = 'protein'
        if is_nucleotide(str(inseq.seq)):
            alphabet = 'nucleotide'
        seq_info[inseq.id] = {
            'alphabet': alphabet,
            'seq': str(inseq.seq)
        }
    return seq_info

def six_frame_translation(dna_sequence):
    frames = []
    for i in range(3):  # forward 3 frames
        seq_mod = dna_sequence[i:]  
        if len(seq_mod) % 3 != 0: 
            seq_mod += "N" * (3 - len(seq_mod) % 3)  # add N to seq
        translated_seq = seq_mod.translate()
        frames.append(translated_seq)
    
    for i in range(3):  # reverse comp 3 frames 
        seq_mod = dna_sequence.reverse_complement()[i:]
        if len(seq_mod) % 3 != 0:
            seq_mod += "N" * (3 - len(seq_mod) % 3)
        translated_seq = seq_mod.translate()
        frames.append(translated_seq)
    
    return frames


# put hmmscan alignment detail to hash 
def aln_to_hash(hmmscan_detail, ga_cutoff):
    aln_hash = {}  # key: protein seq id; value: alignment information
    detail_line = hmmscan_detail.strip().split('\n')
    for line in detail_line:
        a = line.split('\t')
        pfam_id, score, evalue = a[1], float(a[9]), float(a[10])
        pfam_id = re.sub(r'\..*', '', pfam_id)
        if pfam_id in ga_cutoff and score < ga_cutoff[pfam_id]:
            continue
        if pfam_id not in ga_cutoff and evalue > 1e-3:
            continue

        if a[0] in aln_hash:
            aln_hash[a[0]] += line + "\n"
        else:
            aln_hash[a[0]] = line + "\n"
    return aln_hash

# load plantsp family description to hash 
def pk_to_hash(file):
    hash_ = {}
    try:
        with open(file, 'r') as pfh:
            for line in pfh:
                pm = line.strip().split('\t', 2)
                hash_[pm[0]] = pm[1]
    except IOError as e:
        raise IOError(f"Can not open protein kinase description file: {file} {e}")
    return hash_


def split_domain_num(domain_num):
    # Check if the input format is correct, with a '#' character
    if '#' not in domain_num:
        raise ValueError(f"[ERR]domain num format 1 {domain_num}")

    # Split the input string by '#'
    parts = domain_num.split('#')

    # Check if the input string was correctly split into two parts
    if len(parts) != 2:
        raise ValueError(f"[ERR]domain num format 2 {domain_num}")

    # Check if the second part is a positive integer
    if not parts[1].isdigit() or int(parts[1]) <= 0:
        raise ValueError(f"[ERR]domain num format 3 {domain_num}")

    # Repeat the domain ID as many times as specified, separated by commas
    domain_id = ','.join([parts[0]] * int(parts[1]))

    return domain_id


def parse_domain_rule(domain_rule):
    if domain_rule == 'NA':
        return {}

    domain_combination = {}
    r = domain_rule.split(';')

    for rule in r:
        domain_combination_sub = {}

        m = rule.split('--')

        for member in m:
            p = member.split(':')

            if len(p) < 1:
                raise ValueError(f"[ERR] {member}")

            # For the first domain combination sub
            if len(domain_combination_sub) == 0:
                for domain in p:
                    domain_id = split_domain_num(domain)
                    domain_combination_sub[domain_id] = 1
                continue

            # For the single domain in this member
            if len(p) == 1:
                domain_id = split_domain_num(p[0])
                for com in sorted(domain_combination_sub.keys()):
                    del domain_combination_sub[com]
                    com += f",{domain_id}"
                    domain_combination_sub[com] = 1
                continue

            # For the multiply domains in this member
            if len(p) > 1:
                for com in sorted(domain_combination_sub.keys()):
                    del domain_combination_sub[com]
                    for domain in p:
                        domain_id = split_domain_num(domain)
                        new_com = f"{com},{domain_id}"
                        domain_combination_sub[new_com] = 1

        # Put sub domain combination to domain combination
        for com in sorted(domain_combination_sub.keys()):
            domain_combination[com] = 1

    return domain_combination

'''
# Example:
# rules = load_rule('path_to_your_rule_file.txt')
# print(rules)
'''
def load_rule(rule_file):
    rule_obj = {}
    
    with open(rule_file, 'r') as fh:
        id, name, family, required, auxiliary, forbidden, type, desc = '', '', '', '', '', '', '', ''
        
        for line in fh:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            
            if line.startswith('//'):  # New rule starts
                if not all([id, name, family, required, auxiliary, forbidden, type, desc]):
                    raise ValueError(f"Error: Undefined rule member {id}")
                rule_obj[id] = {
                    'name': name,
                    'family': family,
                    'required': parse_domain_rule(required),
                    'auxiliary': parse_domain_rule(auxiliary),
                    'forbidden': parse_domain_rule(forbidden),
                    'type': type,
                    'desc': desc
                }
                id, name, family, required, auxiliary, forbidden, type, desc = '', '', '', '', '', '', '', ''
            elif line.startswith('ID:'):
                id = line.replace('ID:', '').strip()
            elif line.startswith('Name:'):
                name = line.replace('Name:', '').strip()
            elif line.startswith('Family:'):
                family = line.replace('Family:', '').strip()
            elif line.startswith('Required:'):
                required = line.replace('Required:', '').strip()
            elif line.startswith('Auxiliary:'):
                auxiliary = line.replace('Auxiliary:', '').strip()
            elif line.startswith('Forbidden:'):
                forbidden = line.replace('Forbidden:', '').strip()
            elif line.startswith('Type:'):
                type = line.replace('Type:', '').strip()
            elif line.startswith('Desc:'):
                desc = line.replace('Desc:', '').strip()
    
    return rule_obj


def print_rule(rule_pack):
    for id in sorted(rule_pack.keys()):
        print(id)
        print(rule_pack[id]['name'])
        print(rule_pack[id]['family'])
        print(rule_pack[id]['type'])
        print(rule_pack[id]['desc'])

        print("Required:")
        for d in sorted(rule_pack[id]['required']):
            print(d)

        print("Auxiliary:")
        for d in sorted(rule_pack[id]['auxiliary']):
            print(d)

        print("Forbidden:")
        for d in sorted(rule_pack[id]['forbidden']):
            print(d)

    sys.exit()

def load_ga_cutoff(pfam_db, correct_ga, sfam_db):
    ga_cutoff = {}
    pfam_id = ''
    ga_score = float()

    # Load GA cutoff from sfam_db
    with open(sfam_db) as fh0:
        for line in fh0:
            line = line.strip()
            if line.startswith('ACC'):
                pfam_id = line.split()[1].split('.')[0]
            elif line.startswith('GA'):
                ga_score = float(line.split()[2].rstrip(';'))
            elif line == '//':
                if not pfam_id:
                    print("[WARN]no pfam id")
                if not ga_score:
                    print(f"[WARN]no ga score {pfam_id}")
                ga_cutoff[pfam_id] = ga_score
                pfam_id = ''
                ga_score = float()

    # Load GA cutoff from pfam_db
    with open(pfam_db) as fh1:
        for line in fh1:
            line = line.strip()
            if line.startswith('ACC'):
                pfam_id = line.split()[1].split('.')[0]
            elif line.startswith('GA'):
                ga_score = float(line.split()[2].rstrip(';'))
            elif line == '//':
                if not pfam_id:
                    print("[WARN]no pfam id")
                if not ga_score:
                    print(f"[WARN]no ga score {pfam_id}")
                ga_cutoff[pfam_id] = ga_score
                pfam_id = ''
                ga_score = float()

    # Correct GA cutoff
    stop = False
    with open(correct_ga) as fh2:
        for line in fh2:
            line = line.strip()
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                continue
            pfam_id, ga_score = parts
            pfam_id = pfam_id.split('.')[0]
            ga_score = float(ga_score.rstrip(';'))
            if pfam_id not in ga_cutoff:
                print(f"[ERR]no correct pfam id  {line}")
                stop = True
            if not ga_score:
                print(f"[ERR]no correct ga score {line}")
                stop = True
            ga_cutoff[pfam_id] = ga_score

    if stop:
        sys.exit()

    # Uncomment the following line to print the count of GA scores loaded
    # print(f"{len(ga_cutoff)} record has GA score")
    
    return ga_cutoff


# functions for tf classification 
# -- compare two arrays
#    Example usage:
#    array_A = ['pfam001', 'pfam002', 'pfam003']
#    score = [100, 200, 300]
#    array_B = ['pfam001', 'pfam004']
#    print(compare_array(array_A, score, array_B))
def compare_array(array_A, score, array_B):
    a = list(array_A)
    b = list(array_B)
    s = list(score)

    # Convert list to dictionary
    # ua: key: pfam_id, value: num
    # ub: key: required_pfam_id, value: num
    # sa: key: pfam_id, value: score list
    ua = {}
    ub = {}
    sa = {}

    for idx, item in enumerate(a):
        ua[item] = ua.get(item, 0) + 1
        sa[item] = sa.get(item, []) + [s[idx]]

    for item in b:
        ub[item] = ub.get(item, 0) + 1

    # Compare two dictionaries
    # Find the best score for match
    match = 0
    match_score = 0

    for d in sorted(ub.keys()):
        if d in ua and ua[d] >= ub[d]:
            match += 1
            sorted_scores = sorted(sa[d], reverse=True)
            n = ub[d]
            for score in sorted_scores:
                match_score += score
                n -= 1
                if n == 0:
                    break

    # Return match status
    # 0, do not match
    # 1, partially match
    # 2, full match
    match_status = 0
    if match > 0:
        match_status = 1
    if match == len(ub):
        match_status = 2

    return match_status, match_score

# functions for tf classification 
# -- compare rule
def compare_rule(hmm_hit, hmm_score, hmm_hit_s, hmm_score_s, rule_pack):
    rule = rule_pack  # Unpack the rule

    # Compare the hits with rules, including required, auxiliary, and forbidden domains
    # The comparison will return match status:
    # 0, do not match
    # 1, partially match
    # 2, full match
    # The assign family to each protein according to hits and rules
    rule_id = 'NA'

    hits = hmm_hit.split('\t')
    score = list(map(float, hmm_score.split('\t')))
    hits_s = hmm_hit_s.split('\t')
    score_s = list(map(float, hmm_score_s.split('\t')))

    if len(hits) != len(score):
        raise ValueError("[ERR]hit num do not match score num")

    total_domain = 0  # Number of required domains used in classification
    total_score = 0  # Total of score of required domains

    for rid in sorted(rule.keys()):
        required_h = rule[rid]['required']
        auxiliary_h = rule[rid]['auxiliary']
        forbidden_h = rule[rid]['forbidden']

        # Compare forbidden with hits
        f_status = 0
        if forbidden_h != 'NA':
            for forbidden in sorted(forbidden_h.keys()):
                f = forbidden.split(',')
                match_status, match_score = compare_array(hits, score, f)
                if match_status > 0:
                    f_status = 1
        if f_status == 1:
            continue

        # Compare required with hits
        r_status = 0
        for required in sorted(required_h.keys()):
            r = required.split(',')
            match_status, match_score = compare_array(hits, score, r)

            # Define several families need compare in sequence level
            match_status_b, match_score_b = (0, 0)
            if rid in ['T0008', 'T0009', 'T0011', 'T0023']:
                match_status_b, match_score_b = compare_array(hits_s, score_s, r)
            domain_num = len(r)
            if match_status == 2:
                r_status = 1

            if match_status == 2 or match_status_b == 2:
                if rid == 'T9999':
                    if total_domain == 0 and total_score == 0:
                        total_domain = domain_num
                        total_score = match_score
                        rule_id = rid
                    else:
                        if domain_num > total_domain and match_score > total_score:
                            total_domain = domain_num
                            total_score = match_score
                            rule_id = rid
                else:
                    if total_domain == 0 and total_score == 0:
                        total_domain = domain_num
                        total_score = match_score
                        rule_id = rid
                    else:
                        if domain_num >= total_domain and match_score > total_score:
                            total_domain = domain_num
                            total_score = match_score
                            rule_id = rid

        # The auxiliary comparison code is commented out as it may be used in the future

    return rule_id

# functions for tf classification 
# -- tf identification & classification 
def itak_tf_identify(hmmscan_hit, hmmscan_detail, hmmscan_hit_b, ga_cutoff, tf_rule):
    # Remove trailing newline characters
    hmmscan_hit = hmmscan_hit.rstrip()
    hmmscan_detail = hmmscan_detail.rstrip()
    hmmscan_hit_b = hmmscan_hit_b.rstrip()

    # Dictionary to store result with query id as key and tid of TF as value
    qid_tid = {}

    # Dictionary to store all hit domains with query ID and score as key and value
    query_hits_all = {}
    query_hits = {}

    # Parse hmmscan hits
    for line in hmmscan_hit.split('\n'):
        query_id, pfam_id, score, evalue = line.split('\t')
        pfam_id = pfam_id.split('.')[0] if 'PF' in pfam_id else pfam_id
        if pfam_id not in ga_cutoff:
            raise Exception(f"[ERR]undef GA score for {pfam_id}")
        if float(score) < ga_cutoff[pfam_id]:
            continue

        if query_id in query_hits:
            query_hits[query_id]['pid'] += f"\t{pfam_id}"
            query_hits[query_id]['score'] += f"\t{score}"
        else:
            query_hits[query_id] = {'pid': pfam_id, 'score': score}
            query_hits_all[query_id] = 1

    # Dictionary to store hit domains of sequence with query ID and score as key and value
    query_hits_s = {}

    # Parse hmmscan hits for sequences
    for line in hmmscan_hit_b.split('\n'):
        query_id, pfam_id, score, evalue = line.split('\t')
        pfam_id = pfam_id.split('.')[0] if 'PF' in pfam_id else pfam_id
        if pfam_id not in ga_cutoff:
            raise Exception(f"[ERR]undef GA score for {pfam_id}")
        if float(score) < ga_cutoff[pfam_id]:
            continue

        if query_id in query_hits_s:
            query_hits_s[query_id]['pid'] += f"\t{pfam_id}"
            query_hits_s[query_id]['score'] += f"\t{score}"
        else:
            query_hits_s[query_id] = {'pid': pfam_id, 'score': score}
            query_hits_all[query_id] = 1

    # Compare hits with rules
    for qid in sorted(query_hits_all.keys()):
        hits = query_hits[qid]['pid'] if qid in query_hits else ''
        score = query_hits[qid]['score'] if qid in query_hits else ''

        hits_s = query_hits_s[qid]['pid'] if qid in query_hits_s else ''
        score_s = query_hits_s[qid]['score'] if qid in query_hits_s else ''

        if hits and score and hits_s and score_s:
<<<<<<< HEAD
            rule_id = compare_rule(hits, score, hits_s, score_s, tf_rule)
=======
	    rule_id = compare_rule(hits, score, hits_s, score_s, tf_rule)
>>>>>>> 7fd45199d6a61dc77072fea061d750ae0a69b0e3
            if rule_id != 'NA':
                qid_tid[qid] = rule_id

    return qid_tid

# functions for tf classification 
# -- Write out tf result to output file
def itak_tf_write_out(qid_tid, seq_info, hmmscan_detail_1, tf_rule, 
                      output_sequence, output_alignment, output_classification):
    # Put hmmscan_detail to dictionary
    q_detail = {}
    hmmscan_detail_1 = hmmscan_detail_1.strip()
    for a in hmmscan_detail_1.split('\n'):
        b = a.split('\t')
        q_detail.setdefault(b[0], '')
        q_detail[b[0]] += a + "\n"

    # Write to output files
    with open(output_sequence, 'w') as out1, \
         open(output_alignment, 'w') as out2, \
         open(output_classification, 'w') as out3:

        for qid in sorted(qid_tid.keys()):
            tid = qid_tid[qid]
            tname = tf_rule[tid]['name']
            tfamily = tf_rule[tid]['family']
            type_ = tf_rule[tid]['type']
            desc = tf_rule[tid]['desc']
            qseq = seq_info[qid]['seq']
            align = q_detail.get(qid, '')

            out1.write(f">{qid} [{type_}]{tname}:{tfamily}--{desc}\n{qseq}\n")
            out2.write(align)
            out3.write(f"{qid}\t{tname}\t{type_}\t{tfamily}\n")


# function for pk classification 
def itak_pk_classify(hmmscan_detail, pkinase_id, other):
    pk_id = dict(pkinase_id)
    hmmscan_detail = hmmscan_detail.strip()
    hit_lines = hmmscan_detail.split('\n')

    # Put hmmscan hit to dictionaries
    hit, score, best_hit_align = {}, {}, {}
    for line in hit_lines:
        a = line.split('\t')
        if a[0] in hit:
            if float(a[9]) > score[a[0]]:
                hit[a[0]] = a[1]
                score[a[0]] = float(a[9])
                best_hit_align[a[0]] = line + "\n"
        else:
            hit[a[0]] = a[1]
            score[a[0]] = float(a[9])
            best_hit_align[a[0]] = line + "\n"
            if a[0] in pk_id:
                del pk_id[a[0]]

    # Classify unaligned pks to other
    for pid in sorted(pk_id.keys()):
        hit[pid] = other
        best_hit_align[pid] = ''

    return hit, best_hit_align

"""
main function of TF/PK identification.
"""
def itak_identify(args, bin_dir, dbs_dir, database_files):
    # options, files


    # itak_identify(args, bin_path, dbs_path, database_files)
    # check and set up output folder
    for f in args.seq_files:
        output_dir = f"{f}_output"
        if args.output is not None:
            output_dir = args.output
        temp_dir = f"{f}_temp"
        if os.path.exists(output_dir):
            print(f"[WARN]output folder exist: {output_dir}")
        if os.path.exists(temp_dir):
            print(f"[WARN]temp folder exist: {temp_dir}")
        if not os.path.isfile(f) or os.path.getsize(f) == 0:
            print("[ERR]input file not exist")
            sys.exit(1)

    # set number of processes
    cpu = '20'
    if args.process > 0:
        cpu = str(args.process)

    # Frame for translate
    frame = '6'  # Default is 6
    if args.frame in ('3F', '3R'):
        frame = args.frame

    mode = args.mode

    # set vars for database files, rule files, and hmmscan
    tfam_db = dbs_dir + "/Tfam_domain.hmm"  # database for transcription factors (subset of Pfam-A + customized). for quick mode
    sfam_db = dbs_dir + "/TF_selfbuild.hmm"  # self-build -- backup
    plantsp_db = dbs_dir + "/PlantsPHMM3_89.hmm"  # plantsP kinase
    shiu_db = dbs_dir + "/Plant_Pkinase_fam.hmm"  # shiu kinase database
    Psub_wnk1 = dbs_dir + "/Pkinase_sub_WNK1.hmm"  # wnk1 hmm  -- for plantsp 
    Psub_MAK = dbs_dir + "/Pkinase_sub_MAK.hmm"  # MAK hmm  -- for plantsp
    pfam_db = dbs_dir + "/Pfam-A.hmm"  # Pfam-A -- for normal mode, obsolete

    tf_rule = dbs_dir + "/TF_Rule.txt"  # Rules for Transcription Factors
    correct_ga = dbs_dir + "/GA_table.txt"  # update GA cutoff
    pk_desc = dbs_dir + "/PK_class_desc.txt"  # PK family description (for PPC)

    hmmscan_bin = bin_dir + "/hmmscan"  # hmmscan

    # Load rules, GA cutoffs, and PK descriptions into dictionaries
    # -- load tf rules
    tf_rule = load_rule(os.path.join(dbs_dir, "TF_Rule.txt"))
    
    # -- Initialize an empty dictionary for ga_cutoff, then Load ga_cutoff based on the mode
    ga_cutoff = {}
    if mode == 'normal':
        ga_cutoff = load_ga_cutoff(pfam_db, correct_ga, sfam_db)
    elif mode == 'quick':
        ga_cutoff = load_ga_cutoff(tfam_db, correct_ga, sfam_db)
    else:
        print(f"Mode error: {mode}")
        sys.exit(1)

    # -- Convert pk_desc to a hash (dictionary in Python)
    pkid_des = pk_to_hash(pk_desc)

	# +++++ prepare files for norm mode (obsolete) +++++
    """
    if ($mode eq 'normal') {
        my $norm_db_cmd = "# please download and prepare Pfam database using below command:\n".
            "wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz\n".
            "gunzip  Pfam-A.hmm.gz\n".
            "mv Pfam-A.hmm database\n".
            "$hmmpress_bin database/Pfam-A.hm\n\n".
            "# please prepare hmm database for self-build domain\n".
            "$hmmpress_bin $sfam_db\n\n";

        print $norm_db_cmd and exit unless -s $pfam_db;
        die "[ERR]no self-build file $sfam_db\n" unless -s $sfam_db;

        foreach my $db (($pfam_db, $sfam_db)) {
            unless (-s $db.".h3f" && -s $db.".h3i" && -s $db.".h3m" && -s $db.".h3p") {
                warn "[WARN]no database file $db\n";
                run_cmd("$hmmpress_bin -f $db");
            }
        }
    }
    """

    # Main processing loop for input files
    for f in args.seq_files:
        temp_dir = f"{f}_temp"
        output_dir = f"{f}_output"
        if args.output is not None:
            output_dir = args.output
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        input_protein_f = os.path.join(temp_dir, "protein_seq.fa")
        tmp_pfam_hmmscan = os.path.join(temp_dir, "protein_seq.pfam.hmmscan.txt")

        # Initialize report
        report_info = '#' * 80

        # Load sequences from a file into a dictionary
        seq_info = seq_to_hash(f)
        report_info += f"\nLoad {len(seq_info)} sequences from {f}\n"

        protein_num, nucleotide_num = 0, 0

        # Open the output file for writing protein sequences
        with open(input_protein_f, 'w') as outp:
            for id in sorted(seq_info.keys()):

                if seq_info[id]['alphabet'] == 'protein':
                    protein_num += 1
                    outp.write(f">{id}\n{seq_info[id]['seq']}\n")
                else:
                    nucleotide_num += 1
                    # just print nucleotide to screen 
                    print(f">{id}\n{seq_info[id]['seq']}\n")

                    # Translate to proteins  
                    translated_frames = six_frame_translation(Seq(seq_info[id]['seq']))
                    for i, frame in enumerate(translated_frames, start=1):
                        new_id = id + '-' + str(i)
                        outp.write(f">{new_id}\n{frame}\n")
                        seq_info[new_id] = {'seq': frame}
        report_info += f"  {protein_num} of input sequences are protein\n  {nucleotide_num} of input sequences are nucleotide\n"

        # ==== Part A TF identification ====
        # ==== A1. compare input seq with database ====
        hmmscan_hit_1 = hmmscan_detail_1 = hmmscan_hit_1b = None

        if mode == 'normal':
            hmmscan_command1 = f"{hmmscan_bin} --acc --notextw --cpu {cpu} -o {tmp_pfam_hmmscan}.a {pfam_db} {input_protein_f}"
            hmmscan_command2 = f"{hmmscan_bin} --acc --notextw --cpu {cpu} -o {tmp_pfam_hmmscan}.b {sfam_db} {input_protein_f}"
            
            if not os.path.getsize(tmp_pfam_hmmscan + ".a"):
                run_cmd(hmmscan_command1)
            if not os.path.getsize(tmp_pfam_hmmscan + ".b"):
                run_cmd(hmmscan_command2)
            
            hmmscan_hit_1, hmmscan_detail_1, hmmscan_hit_1b = itakm.parse_hmmscan_result(tmp_pfam_hmmscan + ".a")
            hmmscan_hit_2, hmmscan_detail_2, hmmscan_hit_2b = itakm.parse_hmmscan_result(tmp_pfam_hmmscan + ".b")
            
            hmmscan_hit_1 += hmmscan_hit_2
            hmmscan_detail_1 += hmmscan_detail_2
            hmmscan_hit_1b += hmmscan_hit_2b
        else:
            hmmscan_command = f"{hmmscan_bin} --acc --notextw --cpu {cpu} -o {tmp_pfam_hmmscan} {tfam_db} {input_protein_f}"
            # print(hmmscan_command)
            if os.path.exists(tmp_pfam_hmmscan) and os.path.getsize(tmp_pfam_hmmscan) > 0:
                print("alignment file exist")
            else:
                run_cmd(hmmscan_command)
            
            hmmscan_hit_1, hmmscan_detail_1, hmmscan_hit_1b = itakm.parse_hmmscan_result(tmp_pfam_hmmscan)

        align_pfam_hash = aln_to_hash(hmmscan_detail_1, ga_cutoff)

        # Part A2 TF identification
        qid_tid = itak_tf_identify(hmmscan_hit_1, hmmscan_detail_1, hmmscan_hit_1b, ga_cutoff, tf_rule)
        report_info += f"  {len(qid_tid)} of proteins were identified as transcription factors or transcriptional regulators\n"

        # Part A3 save the result
        output_sequence = f"{output_dir}/tf_sequence.fasta"
        output_alignment = f"{output_dir}/tf_alignment.txt"
        output_classification = f"{output_dir}/tf_classification.txt"
        itak_tf_write_out(qid_tid, seq_info, hmmscan_detail_1, tf_rule, output_sequence, output_alignment, output_classification)

        print(report_info)

        # ==== Part X Specific Family Identification ==== 
        spec_families = check_specific_families(dbs_dir)
        if args.specific in spec_families:
            # ==== X1. compare input seq with specific database ====
            # prepare rules and db   
            spec_rule = load_rule(os.path.join(dbs_dir, "specific", args.specific + ".txt"))
            spec_hmm = os.path.join(dbs_dir, "specific", args.specific + ".hmm")
            spec_ga_cutoff = load_ga_cutoff(tfam_db, correct_ga, spec_hmm)

            if not all(os.path.isfile(f"{spec_hmm}.{ext}") for ext in ["h3f", "h3i", "h3m", "h3p"]):
                print(f"[WARN]no database file {spec_hmm}")
                subprocess.call([os.path.join(bin_dir, "hmmpress"), "-f", spec_hmm])
            else:
                print(f"Database file {spec_hmm} is ready.")

            # temp file for spec hmmscan result 
            spec_tmp_hmmscan = os.path.join(temp_dir, args.specific + ".protein_seq.hmmscan.txt")

            spec_hmmscan_command = f"{hmmscan_bin} --acc --notextw --cpu {cpu} -o {spec_tmp_hmmscan} {spec_hmm} {input_protein_f}"
            # print(hmmscan_command)
            if os.path.exists(spec_tmp_hmmscan) and os.path.getsize(spec_tmp_hmmscan) > 0:
                print(f"{args.specific} alignment file exist")
            else:
                run_cmd(spec_hmmscan_command)
            
            hmmscan_hit_3, hmmscan_detail_3, hmmscan_hit_3b = itakm.parse_hmmscan_result(spec_tmp_hmmscan)

            if not hmmscan_hit_3:
                report_info = f"  The input sequences are not recognized and classified by rules: {args.specific} \n"
            else:
                # Part X2 identification
                spec_qid_tid = itak_tf_identify(hmmscan_hit_3, hmmscan_detail_3, hmmscan_hit_3b, spec_ga_cutoff, spec_rule)
                report_info = f"  {len(spec_qid_tid)} of proteins were processed with specific rules: {args.specific} \n"

                # Part X3 save the result
                spec_output_sequence = f"{output_dir}/{args.specific}_sequence.fasta"
                spec_output_alignment = f"{output_dir}/{args.specific}_alignment.txt"
                spec_output_classification = f"{output_dir}/{args.specific}_classification.txt"
                itak_tf_write_out(spec_qid_tid, seq_info, hmmscan_detail_3, spec_rule, spec_output_sequence, spec_output_alignment, spec_output_classification)
            print(report_info)
            
        else:
            if args.specific is not None:
                print(f"Can not find rules and database for {args.specific}")

        # ==== Part B PK identification ====
        # ==== B1. get protein kinase seqs ====
        pkinase_id = {}
        hit_lines = hmmscan_hit_1.strip().split('\n')
        for line in hit_lines:
            a = line.split('\t')
            a[1] = a[1].split('.')[0]  # Equivalent to $a[1] =~ s/\..*//;
            
            if a[1] in ['PF00069', 'PF07714']:
                if a[1] not in ga_cutoff:
                    raise Exception(f"[ERR]no GA for {a[1]}")
                if float(a[2]) >= ga_cutoff[a[1]]:
                    pkinase_id[a[0]] = 1

        if len(pkinase_id) == 0:
            report_info = "  no protein was identified as protein kinase\n"
            report_info += "Finished\n"
            report_info += "#" * 80
            print(report_info)

            if not debug:
                if os.path.getsize(temp_dir) > 0:
                    os.system(f"rm -rf {temp_dir}")
            
            #if 'z' in options:
            #    os.system(f"tar -czf {output_dir}.tgz {output_dir}")
            #if 's' in options:
            #    send_mail(options['s'], f)  # Assuming send_mail is a function you have defined elsewhere
            # not protein kinase seq, parse next input sequence file
            continue
        
        tmp_pkinase_seq = os.path.join(temp_dir, "pkinase_seq.fa")
        with open(tmp_pkinase_seq, 'w') as out1:
            for id in sorted(pkinase_id.keys()):
                out1.write(f">{id}\n{seq_info[id]['seq']}\n")

        # ==== B2. compare protein kinase seqs with database ====
        tmp_plantsp_hmmscan = os.path.join(temp_dir, "protein_seq.plantsp.hmmscan.txt")
        tmp_shiu_hmmscan = os.path.join(temp_dir, "protein_seq.shiu.hmmscan.txt")
        # tmp_rkd_hmmscan = os.path.join(temp_dir, "protein_seq.rkd.hmmscan.txt")

        plantsp_hmmscan_cmd = f"{hmmscan_bin} --acc --notextw --cpu {cpu} -o {tmp_plantsp_hmmscan} {plantsp_db} {tmp_pkinase_seq}"
        shiu_hmmscan_cmd = f"{hmmscan_bin} --acc --notextw --cpu {cpu} -o {tmp_shiu_hmmscan} {shiu_db} {tmp_pkinase_seq}"
        # rkd_hmmscan_cmd = f"{hmmscan_bin} --acc --notextw --cpu {cpu} -o {tmp_rkd_hmmscan} {rkd_db} {tmp_pkinase_seq}"

        
        # === the plantsp classification will run with parameter c ===
        plantsp_hit, plantsp_detail, plantsp_hit_b = '', '', ''
        plantsp_cat, plantsp_aln = '', ''
        if args.classify is True:
            if not os.path.exists(tmp_plantsp_hmmscan):
                run_cmd(plantsp_hmmscan_cmd)
            
            plantsp_hit, plantsp_detail, plantsp_hit_b = itakm.parse_hmmscan_result(tmp_plantsp_hmmscan)
            plantsp_cat, plantsp_aln = itak_pk_classify(plantsp_detail, pkinase_id, "PPC:5.2.1")
            
            # Classification of sub pkinase
            sub_classifications = [
                (os.path.join(dbs_dir, "Pkinase_sub_WNK1.hmm"), "30", "PPC:4.1.5", "PPC:4.1.5.1"),
                (os.path.join(dbs_dir, "Pkinase_sub_MAK.hmm"), "460.15", "PPC:4.5.1", "PPC:4.5.1.1")
            ]
            
            for sub in sub_classifications:
                if len(sub) != 4:
                    raise Exception(f"[ERR]sub classify info {','.join(sub)}")
                
                hmm_profile, cutoff, cat, sub_cat = sub
                
                seq_num = 0
                ppc_seq = os.path.join(temp_dir, "temp_ppc_seq")
                with open(ppc_seq, 'w') as ppcfh:
                    for seq_id in sorted(plantsp_cat.keys()):
                        if plantsp_cat[seq_id] == cat:
                            if seq_id not in seq_info or 'seq' not in seq_info[seq_id]:
                                raise Exception(f"[ERR]seq id: {seq_id}")
                            ppcfh.write(f">{seq_id}\n{seq_info[seq_id]['seq']}\n")
                            seq_num += 1
                
                if seq_num == 0:
                    continue
                print(f"{seq_num}\t{cat}")
                
                ppc_hmm_result = os.path.join(temp_dir, "temp_ppc_sub_hmmscan.txt")
                hmm_cmd = f"{hmmscan_bin} --acc --notextw --cpu {cpu} -o {ppc_hmm_result} {hmm_profile} {ppc_seq}"
                run_cmd(hmm_cmd)
                
                ppc_hits, ppc_detail, ppc_hits_b = itakm.parse_hmmscan_result(ppc_hmm_result)
                for hit in ppc_detail.split('\n'):
                    a = hit.split('\t')
                    if float(a[9]) >= float(cutoff):
                        plantsp_cat[a[0]] = sub_cat
                        plantsp_aln[a[0]] = hit + "\n"

        # === the shiu classification will run by default ===
        if not os.path.exists(tmp_shiu_hmmscan):
            run_cmd(shiu_hmmscan_cmd)

        shiu_hit, shiu_detail, shiu_hit_b = itakm.parse_hmmscan_result(tmp_shiu_hmmscan)

        # ==== B3. PK classification ====	
        shiu_cat, shiu_aln = itak_pk_classify(shiu_detail, pkinase_id, "Group-other")

        # ==== B5 save result =====
        # output plantsp classification
        ppc_cat = os.path.join(output_dir, "PPC_classification.txt")
        ppc_aln = os.path.join(output_dir, "PPC_alignment.txt")

        # option c for report both PPC and Shiu classifcation
        if args.classify is True:
            # print(plantsp_cat)
            with open(ppc_cat, 'w') as ca_fh1, open(ppc_aln, 'w') as al_fh1:
                for pid in sorted(plantsp_cat.keys()):
                    ca_fh1.write(f"{pid}\t{plantsp_cat[pid]}\t{pkid_des[plantsp_cat[pid]]}\n")
                    if pid in align_pfam_hash and pid in plantsp_cat:
                        al_fh1.write(f"{plantsp_aln[pid]}{align_pfam_hash[pid]}")
                    else:
                        raise Exception("Error! Do not have alignments in hmm3 parsed result\n")

        # output Shiu classification
        shiu_cat_file = os.path.join(output_dir, "shiu_classification.txt")
        shiu_aln_file = os.path.join(output_dir, "shiu_alignment.txt")

        with open(shiu_cat_file, 'w') as ca_fh2, open(shiu_aln_file, 'w') as al_fh2:
            for pid in sorted(shiu_cat.keys()):
                ca_fh2.write(f"{pid}\t{shiu_cat[pid]}\n")
                if pid in align_pfam_hash and pid in shiu_cat:
                    al_fh2.write(f"{shiu_aln[pid]}{align_pfam_hash[pid]}")
                else:
                    raise Exception("Error! Do not have alignments in hmm3 parsed result\n")

        # output protein kinase sequences 
        pkinase_seq = os.path.join(output_dir, "pk_sequence.fasta")
        with open(pkinase_seq, 'w') as out_pks:
            for pid in sorted(pkinase_id.keys()):
                cat1 = 'NA'
                cat2 = 'NA'
                if pid in plantsp_cat:
                    cat1 = plantsp_cat[pid]
                if pid in shiu_cat:
                    cat2 = shiu_cat[pid]
                
                if args.classify: # option c for report both PPC and Shiu classifcation
                    out_pks.write(f">{pid} PlantsP:{cat1};Shiu:{cat2}\n{seq_info[pid]['seq']}\n")
                else:
                    out_pks.write(f">{pid} Shiu:{cat2}\n{seq_info[pid]['seq']}\n")

        # report information 
        report_info = f"  {len(pkinase_id)} of proteins were identified as protein kinase\nFinished\n" + "#" * 80
        print(report_info)

        # remove temp folder
        if not debug:
            if os.path.exists(temp_dir) and os.path.isdir(temp_dir):
                os.system(f"rm -rf {temp_dir}")

        # code for online version 
        # if args.compress is True:
        #    os.system(f"tar -czf {output_dir}.tgz {output_dir}")
        # if args.email is not None:
        #    send_mail(args.email, f) 
     
def main():

    # define database files
    database_files = {
        'rule': ["TF_Rule.txt", "GA_table.txt", "PK_class_desc.txt"], 
        'hmm': [
            "Tfam_domain.hmm", "TF_selfbuild.hmm","PlantsPHMM3_89.hmm", 
            "Plant_Pkinase_fam.hmm", "Pkinase_sub_WNK1.hmm", "Pkinase_sub_MAK.hmm"], 
        'other': ["TF_selfbuild.sto"]
        }
    
    # define the url of latest database
    db_url_latest = "https://github.com/kentnf/iTAK/archive/refs/tags/db-v1.tar.gz"

    # check whether the os is linux/macos, whether the hmmer and database files exist
    check_os()
    bin_path = check_hmmer()
    dbs_path = get_db_path()
    
    check_and_prepare_database(dbs_path, db_url_latest, database_files, bin_path)
    # check rule files
    missing_rule = check_db_file_exit(dbs_path, database_files, "rule", bin_path)
    if len(missing_rule) > 0:
        print(f"Can not find rule file: {missing_rule}")
        sys.exit(1)
    
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description='iTAK -- Plant Transcription factor & Protein Kinase Identifier and Classifier')
    parser.add_argument('seq_files', nargs='*' if '-u' in sys.argv or '--update' in sys.argv else '+', help='Input sequence file(s)')
    parser.add_argument('-f', '--frame', help='Translate frame. (3F, 3R, 6; default = 6)', default='6')
    parser.add_argument('-p', '--process', type=int, help='Number of CPUs used for hmmscan. (default = 1)', default=1)
    parser.add_argument('-u', '--update', action='store_true', help='Update the database.')
    parser.add_argument('-s', '--specific', help='User-defined specific family name for identification and classification')
    parser.add_argument('-o', '--output', help='Name of the output directory. (default = \'input file name\' + \'_output\')')
    parser.add_argument('-m', '--mode', help='Mode, quick or normal, please do not change it. (default = quick)', default='quick')
    parser.add_argument('-c', '--classify', action='store_true', help='OBSOLETE: Enable protein kinase PPC classification. (default = disable)')
    #parser.add_argument('-z', '--compress', action='store_true', help='Enable the compression of result files. For online version only. (default = disable)')
    #parser.add_argument('-e', '--email', help='E-mail address. iTAK will send email when analysis is done. For online version only.')
    args = parser.parse_args()
    print(args) # for debug

    # update the database and rules to latest version
    if args.update:
        with tempfile.TemporaryDirectory() as temp_dir:
            # download and extract database files from temp to dbs path
            temp_dir = Path(temp_dir)
            download_and_extract(db_url_latest, temp_dir, dbs_path)

    # identification and classification
    else:
        itak_identify(args, bin_path, dbs_path, database_files)

if __name__ == "__main__":
    main()
