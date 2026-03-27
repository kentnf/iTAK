import os
import re
import sys

from itak.rules import build_domain_hit_collection, compare_rule, rule_attr


def identify_protein_kinases(summary_hits, ga_cutoff):
    pkinase_id = {}
    for hit in summary_hits:
        pfam_id = hit.hit_id.split('.')[0]
        if pfam_id in ['PF00069', 'PF07714']:
            if pfam_id not in ga_cutoff:
                raise Exception(f"[ERR]no GA for {pfam_id}")
            if float(hit.score) >= ga_cutoff[pfam_id]:
                pkinase_id[hit.query_name] = 1
    return pkinase_id


def write_pkinase_sequences(output_dir, pkinase_id, seq_info, plantsp_cat, shiu_cat, include_ppc):
    pkinase_seq = os.path.join(output_dir, "pk_sequence.fasta")
    with open(pkinase_seq, 'w') as out_pks:
        for pid in sorted(pkinase_id.keys()):
            cat1 = plantsp_cat.get(pid, 'NA')
            cat2 = shiu_cat.get(pid, 'NA')
            if include_ppc:
                out_pks.write(f">{pid} PlantsP:{cat1};Shiu:{cat2}\n{seq_info[pid]['seq']}\n")
            else:
                out_pks.write(f">{pid} Shiu:{cat2}\n{seq_info[pid]['seq']}\n")


def aln_to_hash(hmmscan_detail, ga_cutoff):
    aln_hash = {}

    if isinstance(hmmscan_detail, str):
        detail_iter = [line for line in hmmscan_detail.strip().split('\n') if line]
    else:
        detail_iter = hmmscan_detail

    for detail in detail_iter:
        if isinstance(detail, str):
            line = detail
            a = line.split('\t')
            query_id = a[0]
            pfam_id = a[1]
            score = float(a[9])
            evalue = float(a[10])
        else:
            line = detail.to_legacy_line().rstrip('\n')
            query_id = detail.query_name
            pfam_id = detail.hit_id
            score = float(detail.score)
            evalue = float(detail.evalue)

        pfam_id = re.sub(r'\..*', '', pfam_id)
        if pfam_id in ga_cutoff and score < ga_cutoff[pfam_id]:
            continue
        if pfam_id not in ga_cutoff and evalue > 1e-3:
            continue

        if query_id in aln_hash:
            aln_hash[query_id] += line + "\n"
        else:
            aln_hash[query_id] = line + "\n"
    return aln_hash


def load_ga_cutoff(pfam_db, correct_ga, sfam_db):
    ga_cutoff = {}
    pfam_id = ''
    ga_score = float()

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

    return ga_cutoff


def iter_hmmscan_hits(hmmscan_hits):
    if isinstance(hmmscan_hits, str):
        for line in hmmscan_hits.rstrip().split('\n'):
            if not line:
                continue
            query_id, pfam_id, score, evalue = line.split('\t')
            yield query_id, pfam_id, score, evalue
        return

    for hit in hmmscan_hits:
        yield hit.query_name, hit.hit_id, hit.score, hit.evalue


def itak_tf_identify(hmmscan_hit, hmmscan_hit_b, ga_cutoff, tf_rule):
    qid_tid = {}
    query_hits_all = {}
    query_hits = {}

    for query_id, pfam_id, score, evalue in iter_hmmscan_hits(hmmscan_hit):
        pfam_id = pfam_id.split('.')[0] if 'PF' in pfam_id else pfam_id
        if pfam_id not in ga_cutoff:
            raise Exception(f"[ERR]undef GA score for {pfam_id}")
        if float(score) < ga_cutoff[pfam_id]:
            continue

        if query_id in query_hits:
            query_hits[query_id]['pid'].append(pfam_id)
            query_hits[query_id]['score'].append(float(score))
        else:
            query_hits[query_id] = {'pid': [pfam_id], 'score': [float(score)]}
            query_hits_all[query_id] = 1

    query_hits_s = {}
    for query_id, pfam_id, score, evalue in iter_hmmscan_hits(hmmscan_hit_b):
        pfam_id = pfam_id.split('.')[0] if 'PF' in pfam_id else pfam_id
        if pfam_id not in ga_cutoff:
            raise Exception(f"[ERR]undef GA score for {pfam_id}")
        if float(score) < ga_cutoff[pfam_id]:
            continue

        if query_id in query_hits_s:
            query_hits_s[query_id]['pid'].append(pfam_id)
            query_hits_s[query_id]['score'].append(float(score))
        else:
            query_hits_s[query_id] = {'pid': [pfam_id], 'score': [float(score)]}
            query_hits_all[query_id] = 1

    for qid in sorted(query_hits_all.keys()):
        hit_collection = build_domain_hit_collection(
            query_hits[qid]['pid'],
            query_hits[qid]['score'],
        ) if qid in query_hits else build_domain_hit_collection([], [])

        hit_collection_s = build_domain_hit_collection(
            query_hits_s[qid]['pid'],
            query_hits_s[qid]['score'],
        ) if qid in query_hits_s else build_domain_hit_collection([], [])

        if hit_collection.hits and hit_collection_s.hits:
            rule_id = compare_rule(hit_collection, None, hit_collection_s, None, tf_rule)
            if rule_id != 'NA':
                qid_tid[qid] = rule_id

    return qid_tid


def itak_tf_write_out(qid_tid, seq_info, hmmscan_detail_1, tf_rule, output_sequence, output_alignment, output_classification):
    q_detail = {}
    if isinstance(hmmscan_detail_1, str):
        detail_lines = [line for line in hmmscan_detail_1.strip().split('\n') if line]
    else:
        detail_lines = [alignment.to_legacy_line().rstrip('\n') for alignment in hmmscan_detail_1]

    for detail_line in detail_lines:
        columns = detail_line.split('\t')
        q_detail.setdefault(columns[0], '')
        q_detail[columns[0]] += detail_line + "\n"

    with open(output_sequence, 'w') as out1, \
         open(output_alignment, 'w') as out2, \
         open(output_classification, 'w') as out3:
        for qid in sorted(qid_tid.keys()):
            tid = qid_tid[qid]
            tname = rule_attr(tf_rule[tid], 'name')
            tfamily = rule_attr(tf_rule[tid], 'family')
            type_ = rule_attr(tf_rule[tid], 'type')
            desc = rule_attr(tf_rule[tid], 'desc')
            qseq = seq_info[qid]['seq']
            align = q_detail.get(qid, '')

            out1.write(f">{qid} [{type_}]{tname}:{tfamily}--{desc}\n{qseq}\n")
            out2.write(align)
            out3.write(f"{qid}\t{tname}\t{type_}\t{tfamily}\n")


def itak_pk_classify(hmmscan_detail, pkinase_id, other):
    pk_id = dict(pkinase_id)
    if isinstance(hmmscan_detail, str):
        hit_lines = [line for line in hmmscan_detail.strip().split('\n') if line]
    else:
        hit_lines = [alignment.to_legacy_line().rstrip('\n') for alignment in hmmscan_detail]

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

    for pid in sorted(pk_id.keys()):
        hit[pid] = other
        best_hit_align[pid] = ''

    return hit, best_hit_align
