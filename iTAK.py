#!/usr/bin/python3

import sys
import os
import re
import argparse
from Bio import SeqIO
from subprocess import call, CalledProcessError
import smtplib
from email.mime.text import MIMEText
from email.header import Header

import itakm

# Define the iTAK version and debug mode
version = 2.0
debug = False

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
        call(cmd, shell=True)
    except CalledProcessError as e:
        print(f"[ERR]cmd: {cmd}")
        raise

def is_nucleotide(sequence):
    nucleotide_set = set("ACGTUacgtu")
    return set(sequence).issubset(nucleotide_set)

# fasta to hash, key: id, value->alphabet: , value->seq: seq;
def seq_to_hash(input_file):
    """
    put seq info to hash
    """
    from Bio import SeqIO

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

        rule_id = compare_rule(hits, score, hits_s, score_s, tf_rule)
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
def itak_identify(options, files):

    #print(options)
    for f in files:
        output_dir = f"{f}_output"
        if 'o' in options and options['o'] is not None:
            output_dir = options['o']
        temp_dir = f"{f}_temp"
        #print(output_dir)
        #sys.exit()
        if os.path.exists(output_dir):
            print(f"[WARN]output folder exist: {output_dir}")
        if os.path.exists(temp_dir):
            print(f"[WARN]temp folder exist: {temp_dir}")
        if not os.path.isfile(f) or os.path.getsize(f) == 0:
            print("[ERR]input file not exist")
            sys.exit(1)

    # Set CPUs for hmmscan
    cpu = '20'
    if 'p' in options and options['p'] > 0:
        cpu = str(options['p'])

    # Frame for translate
    frame = '6'  # Default is 6
    if 'f' in options and options['f'] in ('3F', '3R'):
        frame = options['f']
    
    mode = options['m']

    # Set database and script paths
    bin_dir = os.path.join(os.getcwd(), "bin")
    dbs_dir = os.path.join(os.getcwd(), "database")
    bin_files = [
        "hmmscan", "hmmpress"
    ]

    rule_files = [
        "TF_Rule.txt", "GA_table.txt", "PK_class_desc.txt"
    ]

    # "Pfam-A.hmm",
    databases = [
        "Tfam_domain.hmm", "TF_selfbuild.hmm",
        "PlantsPHMM3_89.hmm", "Plant_Pkinase_fam.hmm",
        "Pkinase_sub_WNK1.hmm", "Pkinase_sub_MAK.hmm"   
    ]

    tfam_db = dbs_dir + "/Tfam_domain.hmm"  # database for transcription factors (subset of Pfam-A + customized). for quick mode
    pfam_db = dbs_dir + "/Pfam-A.hmm"  # Pfam-A
    sfam_db = dbs_dir + "/TF_selfbuild.hmm"  # self-build
    plantsp_db = dbs_dir + "/PlantsPHMM3_89.hmm"  # plantsP kinase
    shiu_db = dbs_dir + "/Plant_Pkinase_fam.hmm"  # shiu kinase database
    Psub_wnk1 = dbs_dir + "/Pkinase_sub_WNK1.hmm"  # wnk1 hmm
    Psub_MAK = dbs_dir + "/Pkinase_sub_MAK.hmm"  # MAK hmm

    tf_rule = dbs_dir + "/TF_Rule.txt"  # Rules for Transcription Factors
    correct_ga = dbs_dir + "/GA_table.txt"  # update GA cutoff
    pk_desc = dbs_dir + "/PK_class_desc.txt"  # PK family description (for PPC)
    hmmscan_bin = bin_dir + "/hmmscan"  # hmmscan
    hmmpress_bin = bin_dir + "/hmmpress"  # hmmpress

    # Check bin fiels, rule files, and databases
    for f in bin_files:
        file_path = os.path.join(bin_dir, f)
        if not os.path.isfile(file_path) or os.path.getsize(file_path) == 0:
            print(f"[ERR]file not exist: {file_path}")
            sys.exit(1)

    for f in rule_files:
        file_path = os.path.join(dbs_dir, f)
        if not os.path.isfile(file_path) or os.path.getsize(file_path) == 0:
            print(f"[ERR]file not exist: {file_path}")
            sys.exit(1)

    for db in databases:
        db_path = os.path.join(dbs_dir, db)
        if not all(os.path.isfile(f"{db_path}.{ext}") for ext in ["h3f", "h3i", "h3m", "h3p"]):
            print(f"[WARN]no database file {db_path}")
            call([os.path.join(bin_dir, "hmmpress"), "-f", db_path])

    # Load rules, GA cutoffs, and PK descriptions into dictionaries
    # -- load tf rules
    tf_rule = load_rule(os.path.join(dbs_dir, "TF_Rule.txt"))
    
    # -- Initialize an empty dictionary for ga_cutoff, then Load ga_cutoff based on the mode
    ga_cutoff = {}
    if mode == 'normal':
        ga_cutoff = load_ga_cutoff(pfam_db, correct_ga, sfam_db)
    elif mode == 'quick':
        ga_cutoff = load_ga_cutoff(tfam_db, correct_ga, sfam_db)

    # -- Convert pk_desc to a hash (dictionary in Python)
    pkid_des = pk_to_hash(pk_desc)


	# +++++ prepare files for norm mode +++++
	# the normal mode only used for iTAK database
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
    for f in files:
        temp_dir = f"{f}_temp"
        output_dir = f"{f}_output"
        if 'o' in options and options['o'] is not None:
            output_dir = options['o']
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
                    #seqobj = Bio.Seq(seq=seq_info[id]['seq'], id=id)
                    #prots = Bio.SeqUtils.translate_6frames(seqobj)
                    #for i, prot in enumerate(prots):
                    #     nid = prot.id[-3:]
                    #     if nid in frame:
                    #         outp.write(f">{prot.id}\n{prot.seq}\n")
                    #         seq_info[prot.id]['seq'] = prot.seq
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
            
            if 'z' in options:
                os.system(f"tar -czf {output_dir}.tgz {output_dir}")
            if 's' in options:
                send_mail(options['s'], f)  # Assuming send_mail is a function you have defined elsewhere
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
        if 'c' in options and options['c'] is True:
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
        if 'c' in options and options['c'] is True:
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
                
                if 'c' in options: # option c for report both PPC and Shiu classifcation
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
        if 'z' in options and options['z'] is True:
            os.system(f"tar -czf {output_dir}.tgz {output_dir}")
        if 's' in options and options['s'] is not None:
            send_mail(options['s'], f) 

     
def main():

    # Command-line argument parsing
    parser = argparse.ArgumentParser(description='iTAK -- Plant Transcription factor & Protein Kinase Identifier and Classifier')
    parser.add_argument('seq_files', nargs='+', help='Input sequence file(s)')
    parser.add_argument('-f', help='Translate frame. (3F, 3R, 6; default = 6)', default='6')
    parser.add_argument('-p', type=int, help='Number of CPUs used for hmmscan. (default = 1)', default=1)
    parser.add_argument('-o', help='Name of the output directory. (default = \'input file name\' + \'_output\')')
    parser.add_argument('-m', help='Mode, quick or normal, please do not change it. (default = quick)', default='quick')
    parser.add_argument('-c', action='store_true', help='OBSOLETE: Enable protein kinase PPC classification . (default = disable)')
    parser.add_argument('-z', action='store_true', help='Enable the compression of result files. For online version only. (default = disable)')
    parser.add_argument('-s', help='E-mail address. iTAK will send email when analysis is done. For online version only.')
    args = parser.parse_args()
    
    # Convert arguments to dictionary format and call the main function
    options_dict = vars(args)
    print(options_dict)
    itak_identify(options_dict, args.seq_files)

if __name__ == "__main__":
    main()
