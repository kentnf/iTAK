#!/usr/bin/env python

from itak_cli import (
    DB_URL_LATEST,
    DATABASE_FILES,
    db_version,
    debug,
    main,
    run_cli,
    send_mail,
    version,
)
from itak_pipeline import (
    aln_to_hash,
    identify_protein_kinases,
    is_nucleotide,
    iter_hmmscan_hits,
    itak_identify,
    itak_pk_classify,
    itak_tf_identify,
    itak_tf_write_out,
    load_ga_cutoff,
    prepare_input_sequences,
    run_pk_pipeline,
    run_primary_tf_scan,
    run_specific_family_pipeline,
    seq_to_hash,
    six_frame_translation,
    translate_frames,
    translate_three_frames,
    write_pkinase_sequences,
    write_tf_outputs,
)
from itak_rules import (
    build_domain_hit_collection,
    compare_array,
    compare_rule,
    load_rule_records,
    rule_attr,
)
from itak_runtime import (
    AnalysisRequest,
    RuntimeConfig,
    SamplePaths,
    build_analysis_request,
    build_arg_parser,
    build_hmmscan_cmd,
    build_runtime_config,
    build_sample_paths,
    check_db_file_exit,
    check_specific_families,
    download_and_extract,
    ensure_sample_dirs,
    file_has_content,
    get_db_path,
    prepare_runtime_environment,
    resolve_output_dir,
    run_cmd,
    validate_input_files,
)


if __name__ == "__main__":
    main()
