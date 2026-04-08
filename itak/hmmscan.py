import re
from dataclasses import dataclass, field
from typing import List


@dataclass(frozen=True)
class HmmscanSummaryHit:
    query_name: str
    hit_id: str
    score: str
    evalue: str

    def to_legacy_line(self):
        return f"{self.query_name}\t{self.hit_id}\t{self.score}\t{self.evalue}\n"


@dataclass(frozen=True)
class HmmscanAlignment:
    query_name: str
    hit_id: str
    query_start: str
    query_end: str
    hit_start: str
    hit_end: str
    query_seq: str
    alignment: str
    hit_seq: str
    score: str
    evalue: str
    description: str
    query_length: str

    def to_summary_hit(self):
        return HmmscanSummaryHit(self.query_name, self.hit_id, self.score, self.evalue)

    def to_legacy_line(self):
        return (
            f"{self.query_name}\t{self.hit_id}\t{self.query_start}\t{self.query_end}\t"
            f"{self.hit_start}\t{self.hit_end}\t{self.query_seq}\t{self.alignment}\t"
            f"{self.hit_seq}\t{self.score}\t{self.evalue}\t{self.description}\t{self.query_length}\n"
        )


@dataclass(frozen=True)
class HmmscanBestDomainHit:
    query_name: str
    hit_id: str
    score: str
    evalue: str

    def to_legacy_line(self):
        return f"{self.query_name}\t{self.hit_id}\t{self.score}\t{self.evalue}\n"


@dataclass
class HmmscanParseResult:
    summary_hits: List[HmmscanSummaryHit] = field(default_factory=list)
    alignments: List[HmmscanAlignment] = field(default_factory=list)
    best_domain_hits: List[HmmscanBestDomainHit] = field(default_factory=list)

    def to_legacy_tuple(self):
        summary_text = "".join(hit.to_legacy_line() for hit in self.summary_hits)
        alignment_text = "".join(alignment.to_legacy_line() for alignment in self.alignments)
        best_domain_text = "".join(hit.to_legacy_line() for hit in self.best_domain_hits)
        return summary_text, alignment_text, best_domain_text


def merge_parse_results(*results):
    merged = HmmscanParseResult()
    for result in results:
        merged.summary_hits.extend(result.summary_hits)
        merged.alignments.extend(result.alignments)
        merged.best_domain_hits.extend(result.best_domain_hits)
    return merged


def parse_hmmscan_records(hmm_result):
    parse_result = HmmscanParseResult()
    one_result = ""

    try:
        with open(hmm_result, 'r') as rfh:
            for line in rfh:
                if line.startswith("#"):
                    continue
                if line.startswith("//"):
                    _parse_result_block(one_result, parse_result)
                    one_result = ""
                else:
                    one_result += line
    except IOError as e:
        print(f"Cannot open hmmscan result file: {hmm_result} {e}")

    return parse_result


def _parse_result_block(one_result, parse_result):
    if not one_result.strip():
        return

    hits_content = one_result.split(">>")
    hit_head = hits_content[0].split("\n")
    query_name, query_length, no_hits, best_domain_hits = _parse_hit_head(hit_head)
    parse_result.best_domain_hits.extend(best_domain_hits)

    if no_hits:
        return

    for hit_content in hits_content[1:]:
        one_hit = ">>" + hit_content
        alignments = parse_align_records(one_hit, query_name, query_length)
        parse_result.alignments.extend(alignments)
        parse_result.summary_hits.extend(alignment.to_summary_hit() for alignment in alignments)


def _parse_hit_head(hit_head_lines):
    query_name = None
    query_length = None
    no_hits = False
    best_domain_hits = []

    for hit_head_line in hit_head_lines:
        match_query_m = re.search(r"^Query:\s+(\S+)\s+\[M\=(\d+)\]", hit_head_line)
        match_query_l = re.search(r"^Query:\s+(\S+)\s+\[L\=(\d+)\]", hit_head_line)
        match_no_hits = re.search(r"No hits detected that satisfy reporting thresholds", hit_head_line)

        if match_query_m:
            query_name, query_length = match_query_m.groups()
            continue
        if match_query_l:
            query_name, query_length = match_query_l.groups()
            continue
        if match_no_hits:
            no_hits = True
            continue

        parts = hit_head_line.split(None, 11)
        if len(parts) > 8 and re.match(r"^\d+\.\d+", parts[1]):
            best_domain_hits.append(
                HmmscanBestDomainHit(
                    query_name=query_name,
                    hit_id=parts[8],
                    score=parts[1],
                    evalue=parts[0],
                )
            )

    return query_name, query_length, no_hits, best_domain_hits


def parse_align_records(hsp_info, query_name, query_length):
    hsp_lines = hsp_info.strip().split('\n')
    info_by_domain = {}
    hit_id = ""
    hit_desc = ""

    for line in hsp_lines:
        if line.startswith('>> '):
            parts = line.split(None, 2)
            hit_id = parts[1]
            hit_desc = parts[2]
        elif re.match(r'^\s+\d+\s+\W\s+', line):
            info_by_domain[int(line.split()[0])] = line

    alignments = []
    domains = hsp_info.split('== domain')
    for domain_index in range(1, len(domains)):
        domain_lines = domains[domain_index].strip().split('\n')
        info = info_by_domain[domain_index].split()

        query_string = ""
        hit_string = ""
        match_string = ""
        hsp_length = None
        align_pos = None
        match_start = False

        for line in domain_lines:
            if re.match(rf'^\s+{re.escape(hit_id)}\s+\d+\s+\S+\s+\d+', line):
                parts = line.split()
                hit_string = parts[3]
                hsp_length = len(hit_string)
                align_pos = line.find(hit_string)
                match_start = True
                if align_pos < 0:
                    raise ValueError(f"Error! Align start position is negative: {align_pos}\n{line}\n{hit_string}\n")
            elif match_start:
                match_string = line[align_pos:align_pos + hsp_length]
                match_start = False
            elif re.match(rf'^\s+{re.escape(query_name)}', line):
                parts = line.split()
                query_string = parts[3]

        alignments.append(
            HmmscanAlignment(
                query_name=query_name,
                hit_id=hit_id,
                query_start=info[9],
                query_end=info[10],
                hit_start=info[6],
                hit_end=info[7],
                query_seq=query_string,
                alignment=match_string,
                hit_seq=hit_string,
                score=info[2],
                evalue=info[5],
                description=hit_desc,
                query_length=query_length,
            )
        )

    return alignments


def parse_align(hsp_info, query_name, query_length):
    alignments = parse_align_records(hsp_info, query_name, query_length)
    summary_text = "".join(alignment.to_summary_hit().to_legacy_line() for alignment in alignments)
    alignment_text = "".join(alignment.to_legacy_line() for alignment in alignments)
    return summary_text, alignment_text


def parse_format_result(in_file, ga_score_hash):
    out_hash = {}
    hsp_hit_id = {}
    hsp_detail = {}

    in_file_lines = in_file.strip().split('\n')
    uid = 0
    len_padding = len(str(len(in_file_lines)))

    for line in in_file_lines:
        fmm = line.strip().split('\t')
        domain_id = fmm[1].split('.')[0]

        valued = False
        ga_score = ga_score_hash.get(domain_id)

        if ga_score is not None:
            if float(fmm[9]) >= ga_score:
                valued = True
        else:
            if float(fmm[10]) <= 1e-3:
                valued = True

        if valued:
            uid += 1
            zero_padding = str(uid).zfill(len_padding)

            hsp_detail[zero_padding] = line
            hsp_hit_id[zero_padding] = domain_id

            if fmm[0] in out_hash:
                out_hash[fmm[0]] += "\t" + zero_padding
            else:
                out_hash[fmm[0]] = zero_padding

    return out_hash, hsp_hit_id, hsp_detail
