import os
import shutil
import subprocess

from itak.classify import (
    aln_to_hash,
    identify_protein_kinases,
    itak_pk_classify,
    itak_tf_identify,
    itak_tf_write_out,
    load_ga_cutoff,
    write_pkinase_sequences,
)
from itak.hmmscan import merge_parse_results, parse_hmmscan_records
from itak.rules import load_rule_records
from itak.runtime import (
    build_analysis_request,
    build_hmmscan_cmd,
    build_runtime_config,
    build_sample_paths,
    check_specific_families,
    ensure_sample_dirs,
    file_has_content,
    run_cmd,
    validate_input_files,
)
from itak.sequences import prepare_input_sequences


def run_primary_tf_scan(input_protein_f, tmp_pfam_hmmscan, runtime_config):
    if runtime_config.mode == 'normal':
        hmmscan_command1 = build_hmmscan_cmd(
            runtime_config.hmmscan_bin,
            runtime_config.cpu,
            tmp_pfam_hmmscan + ".a",
            runtime_config.pfam_db,
            input_protein_f,
        )
        hmmscan_command2 = build_hmmscan_cmd(
            runtime_config.hmmscan_bin,
            runtime_config.cpu,
            tmp_pfam_hmmscan + ".b",
            runtime_config.sfam_db,
            input_protein_f,
        )

        if not file_has_content(tmp_pfam_hmmscan + ".a"):
            run_cmd(hmmscan_command1)
        if not file_has_content(tmp_pfam_hmmscan + ".b"):
            run_cmd(hmmscan_command2)

        parse_result_1 = parse_hmmscan_records(tmp_pfam_hmmscan + ".a")
        parse_result_2 = parse_hmmscan_records(tmp_pfam_hmmscan + ".b")
        return merge_parse_results(parse_result_1, parse_result_2)

    hmmscan_command = build_hmmscan_cmd(
        runtime_config.hmmscan_bin,
        runtime_config.cpu,
        tmp_pfam_hmmscan,
        runtime_config.tfam_db,
        input_protein_f,
    )
    if file_has_content(tmp_pfam_hmmscan):
        print("alignment file exist")
    else:
        run_cmd(hmmscan_command)

    return parse_hmmscan_records(tmp_pfam_hmmscan)


def write_tf_outputs(output_dir, qid_tid, seq_info, hmmscan_detail, tf_rule):
    itak_tf_write_out(
        qid_tid,
        seq_info,
        hmmscan_detail,
        tf_rule,
        f"{output_dir}/tf_sequence.fasta",
        f"{output_dir}/tf_alignment.txt",
        f"{output_dir}/tf_classification.txt",
    )


def run_specific_family_pipeline(request, paths, seq_info, runtime_config):
    spec_families = check_specific_families(runtime_config.dbs_dir)
    if request.specific not in spec_families:
        if request.specific is not None:
            print(f"Can not find rules and database for {request.specific}")
        return

    spec_rule = load_rule_records(os.path.join(runtime_config.dbs_dir, "specific", request.specific + ".txt"))
    spec_hmm = os.path.join(runtime_config.dbs_dir, "specific", request.specific + ".hmm")
    spec_ga_cutoff = load_ga_cutoff(runtime_config.tfam_db, runtime_config.correct_ga, spec_hmm)

    if not all(os.path.isfile(f"{spec_hmm}.{ext}") for ext in ["h3f", "h3i", "h3m", "h3p"]):
        print(f"[WARN]no database file {spec_hmm}")
        subprocess.run([os.path.join(runtime_config.bin_dir, "hmmpress"), "-f", spec_hmm], check=True)
    else:
        print(f"Database file {spec_hmm} is ready.")

    spec_tmp_hmmscan = os.path.join(paths.temp_dir, request.specific + ".protein_seq.hmmscan.txt")
    spec_hmmscan_command = build_hmmscan_cmd(
        runtime_config.hmmscan_bin,
        runtime_config.cpu,
        spec_tmp_hmmscan,
        spec_hmm,
        paths.input_protein_f,
    )
    if file_has_content(spec_tmp_hmmscan):
        print(f"{request.specific} alignment file exist")
    else:
        run_cmd(spec_hmmscan_command)

    spec_parse_result = parse_hmmscan_records(spec_tmp_hmmscan)
    if not spec_parse_result.summary_hits:
        print(f"  The input sequences are not recognized and classified by rules: {request.specific} \n")
        return

    spec_qid_tid = itak_tf_identify(
        spec_parse_result.summary_hits,
        spec_parse_result.best_domain_hits,
        spec_ga_cutoff,
        spec_rule,
    )
    itak_tf_write_out(
        spec_qid_tid,
        seq_info,
        spec_parse_result.alignments,
        spec_rule,
        f"{paths.output_dir}/{request.specific}_sequence.fasta",
        f"{paths.output_dir}/{request.specific}_alignment.txt",
        f"{paths.output_dir}/{request.specific}_classification.txt",
    )
    print(f"  {len(spec_qid_tid)} of proteins were processed with specific rules: {request.specific} \n")


def run_pk_pipeline(request, paths, seq_info, primary_parse_result, align_pfam_hash, runtime_config):
    pkinase_id = identify_protein_kinases(primary_parse_result.summary_hits, runtime_config.ga_cutoff)
    if len(pkinase_id) == 0:
        return "  no protein was identified as protein kinase\nFinished\n" + "#" * 80

    tmp_pkinase_seq = os.path.join(paths.temp_dir, "pkinase_seq.fa")
    with open(tmp_pkinase_seq, 'w') as out1:
        for seq_id in sorted(pkinase_id.keys()):
            out1.write(f">{seq_id}\n{seq_info[seq_id]['seq']}\n")

    tmp_plantsp_hmmscan = os.path.join(paths.temp_dir, "protein_seq.plantsp.hmmscan.txt")
    tmp_shiu_hmmscan = os.path.join(paths.temp_dir, "protein_seq.shiu.hmmscan.txt")
    plantsp_hmmscan_cmd = build_hmmscan_cmd(
        runtime_config.hmmscan_bin,
        runtime_config.cpu,
        tmp_plantsp_hmmscan,
        runtime_config.plantsp_db,
        tmp_pkinase_seq,
    )
    shiu_hmmscan_cmd = build_hmmscan_cmd(
        runtime_config.hmmscan_bin,
        runtime_config.cpu,
        tmp_shiu_hmmscan,
        runtime_config.shiu_db,
        tmp_pkinase_seq,
    )

    plantsp_cat, plantsp_aln = {}, {}
    if request.classify is True:
        if not os.path.exists(tmp_plantsp_hmmscan):
            run_cmd(plantsp_hmmscan_cmd)

        plantsp_parse_result = parse_hmmscan_records(tmp_plantsp_hmmscan)
        plantsp_cat, plantsp_aln = itak_pk_classify(plantsp_parse_result.alignments, pkinase_id, "PPC:5.2.1")

        sub_classifications = [
            (os.path.join(runtime_config.dbs_dir, "Pkinase_sub_WNK1.hmm"), "30", "PPC:4.1.5", "PPC:4.1.5.1"),
            (os.path.join(runtime_config.dbs_dir, "Pkinase_sub_MAK.hmm"), "460.15", "PPC:4.5.1", "PPC:4.5.1.1"),
        ]

        for hmm_profile, cutoff, cat, sub_cat in sub_classifications:
            seq_num = 0
            ppc_seq = os.path.join(paths.temp_dir, "temp_ppc_seq")
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

            ppc_hmm_result = os.path.join(paths.temp_dir, "temp_ppc_sub_hmmscan.txt")
            run_cmd(build_hmmscan_cmd(runtime_config.hmmscan_bin, runtime_config.cpu, ppc_hmm_result, hmm_profile, ppc_seq))

            ppc_parse_result = parse_hmmscan_records(ppc_hmm_result)
            for alignment in ppc_parse_result.alignments:
                if float(alignment.score) >= float(cutoff):
                    plantsp_cat[alignment.query_name] = sub_cat
                    plantsp_aln[alignment.query_name] = alignment.to_legacy_line()

    if not os.path.exists(tmp_shiu_hmmscan):
        run_cmd(shiu_hmmscan_cmd)

    shiu_parse_result = parse_hmmscan_records(tmp_shiu_hmmscan)
    shiu_cat, shiu_aln = itak_pk_classify(shiu_parse_result.alignments, pkinase_id, "Group-other")

    if request.classify is True:
        with open(os.path.join(paths.output_dir, "PPC_classification.txt"), 'w') as ca_fh1, \
             open(os.path.join(paths.output_dir, "PPC_alignment.txt"), 'w') as al_fh1:
            for pid in sorted(plantsp_cat.keys()):
                ca_fh1.write(f"{pid}\t{plantsp_cat[pid]}\t{runtime_config.pkid_des[plantsp_cat[pid]]}\n")
                if pid in align_pfam_hash and pid in plantsp_cat:
                    al_fh1.write(f"{plantsp_aln[pid]}{align_pfam_hash[pid]}")
                else:
                    raise Exception("Error! Do not have alignments in hmm3 parsed result\n")

    with open(os.path.join(paths.output_dir, "shiu_classification.txt"), 'w') as ca_fh2, \
         open(os.path.join(paths.output_dir, "shiu_alignment.txt"), 'w') as al_fh2:
        for pid in sorted(shiu_cat.keys()):
            ca_fh2.write(f"{pid}\t{shiu_cat[pid]}\n")
            if pid in align_pfam_hash and pid in shiu_cat:
                al_fh2.write(f"{shiu_aln[pid]}{align_pfam_hash[pid]}")
            else:
                raise Exception("Error! Do not have alignments in hmm3 parsed result\n")

    write_pkinase_sequences(paths.output_dir, pkinase_id, seq_info, plantsp_cat, shiu_cat, request.classify)
    return f"  {len(pkinase_id)} of proteins were identified as protein kinase\nFinished\n" + "#" * 80


def itak_identify(args, bin_dir, dbs_dir, debug=False):
    request = build_analysis_request(args)
    multiple_inputs = len(request.seq_files) > 1
    validate_input_files(request.seq_files, request.output)
    runtime_config = build_runtime_config(request, bin_dir, dbs_dir, load_ga_cutoff)

    for input_file in request.seq_files:
        paths = build_sample_paths(input_file, request.output, multiple_inputs)
        ensure_sample_dirs(paths)

        report_info = '#' * 80
        seq_info, input_count, protein_num, nucleotide_num = prepare_input_sequences(
            paths.input_file,
            paths.input_protein_f,
            runtime_config.frame,
        )
        report_info += f"\nLoad {input_count} sequences from {paths.input_file}\n"
        report_info += f"  {protein_num} of input sequences are protein\n  {nucleotide_num} of input sequences are nucleotide\n"

        primary_parse_result = run_primary_tf_scan(
            paths.input_protein_f,
            paths.tmp_pfam_hmmscan,
            runtime_config,
        )
        align_pfam_hash = aln_to_hash(primary_parse_result.alignments, runtime_config.ga_cutoff)

        qid_tid = itak_tf_identify(
            primary_parse_result.summary_hits,
            primary_parse_result.best_domain_hits,
            runtime_config.ga_cutoff,
            runtime_config.tf_rule,
        )
        report_info += f"  {len(qid_tid)} of proteins were identified as transcription factors or transcriptional regulators\n"
        write_tf_outputs(paths.output_dir, qid_tid, seq_info, primary_parse_result.alignments, runtime_config.tf_rule)
        print(report_info)

        run_specific_family_pipeline(request, paths, seq_info, runtime_config)

        pk_report = run_pk_pipeline(request, paths, seq_info, primary_parse_result, align_pfam_hash, runtime_config)
        print(pk_report)

        if not debug:
            shutil.rmtree(paths.temp_dir, ignore_errors=True)
