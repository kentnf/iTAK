from importlib import import_module


_EXPORTS = {
    "DB_URL_LATEST": "itak.cli",
    "DATABASE_FILES": "itak.cli",
    "db_version": "itak.cli",
    "debug": "itak.cli",
    "send_mail": "itak.cli",
    "version": "itak.cli",
    "HmmscanAlignment": "itak.hmmscan",
    "HmmscanBestDomainHit": "itak.hmmscan",
    "HmmscanParseResult": "itak.hmmscan",
    "HmmscanSummaryHit": "itak.hmmscan",
    "aln_to_hash": "itak.classify",
    "identify_protein_kinases": "itak.classify",
    "iter_hmmscan_hits": "itak.classify",
    "itak_pk_classify": "itak.classify",
    "itak_tf_identify": "itak.classify",
    "itak_tf_write_out": "itak.classify",
    "load_ga_cutoff": "itak.classify",
    "write_pkinase_sequences": "itak.classify",
    "is_nucleotide": "itak.sequences",
    "prepare_input_sequences": "itak.sequences",
    "seq_to_hash": "itak.sequences",
    "six_frame_translation": "itak.sequences",
    "translate_frames": "itak.sequences",
    "translate_three_frames": "itak.sequences",
    "itak_identify": "itak.pipeline",
    "run_pk_pipeline": "itak.pipeline",
    "run_primary_tf_scan": "itak.pipeline",
    "run_specific_family_pipeline": "itak.pipeline",
    "write_tf_outputs": "itak.pipeline",
    "merge_parse_results": "itak.hmmscan",
    "parse_align": "itak.hmmscan",
    "parse_align_records": "itak.hmmscan",
    "parse_format_result": "itak.hmmscan",
    "parse_hmmscan_records": "itak.hmmscan",
    "DomainHitCollection": "itak.rules",
    "RuleDefinition": "itak.rules",
    "build_domain_hit_collection": "itak.rules",
    "compare_array": "itak.rules",
    "compare_rule": "itak.rules",
    "load_rule_records": "itak.rules",
    "parse_domain_rule": "itak.rules",
    "rule_attr": "itak.rules",
    "split_domain_num": "itak.rules",
    "AnalysisRequest": "itak.runtime",
    "RuntimeConfig": "itak.runtime",
    "SamplePaths": "itak.runtime",
    "build_analysis_request": "itak.runtime",
    "build_arg_parser": "itak.runtime",
    "build_hmmscan_cmd": "itak.runtime",
    "build_runtime_config": "itak.runtime",
    "build_sample_paths": "itak.runtime",
    "check_and_prepare_database": "itak.runtime",
    "check_db_file_exit": "itak.runtime",
    "check_hmmer": "itak.runtime",
    "check_os": "itak.runtime",
    "check_specific_families": "itak.runtime",
    "download_and_extract": "itak.runtime",
    "ensure_sample_dirs": "itak.runtime",
    "file_has_content": "itak.runtime",
    "get_db_path": "itak.runtime",
    "pk_to_hash": "itak.runtime",
    "prepare_runtime_environment": "itak.runtime",
    "resolve_cpu": "itak.runtime",
    "resolve_frame_mode": "itak.runtime",
    "resolve_output_dir": "itak.runtime",
    "run_cmd": "itak.runtime",
    "validate_input_files": "itak.runtime",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module 'itak' has no attribute {name!r}")

    module = import_module(_EXPORTS[name])
    value = getattr(module, name)
    globals()[name] = value
    return value
