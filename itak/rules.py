from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class DomainHitCollection:
    hits: Tuple[str, ...]
    scores: Tuple[float, ...]
    counts: dict
    scores_by_domain: dict


@dataclass(frozen=True)
class RuleDefinition:
    id: str
    name: str
    family: str
    required: dict
    auxiliary: dict
    forbidden: dict
    type: str
    desc: str

    def to_legacy_dict(self):
        return {
            'name': self.name,
            'family': self.family,
            'required': self.required,
            'auxiliary': self.auxiliary,
            'forbidden': self.forbidden,
            'type': self.type,
            'desc': self.desc,
        }


def rule_attr(rule_definition, attr_name):
    if isinstance(rule_definition, RuleDefinition):
        return getattr(rule_definition, attr_name)
    return rule_definition[attr_name]


def split_domain_num(domain_num):
    if '#' not in domain_num:
        raise ValueError(f"[ERR]domain num format 1 {domain_num}")

    parts = domain_num.split('#')
    if len(parts) != 2:
        raise ValueError(f"[ERR]domain num format 2 {domain_num}")

    if not parts[1].isdigit() or int(parts[1]) <= 0:
        raise ValueError(f"[ERR]domain num format 3 {domain_num}")

    return ','.join([parts[0]] * int(parts[1]))


def parse_domain_rule(domain_rule):
    if domain_rule == 'NA':
        return {}

    domain_combination = {}
    for rule in domain_rule.split(';'):
        domain_combination_sub = {}
        for member in rule.split('--'):
            domains = member.split(':')

            if len(domain_combination_sub) == 0:
                for domain in domains:
                    domain_combination_sub[split_domain_num(domain)] = 1
                continue

            if len(domains) == 1:
                domain_id = split_domain_num(domains[0])
                for combination in sorted(domain_combination_sub.keys()):
                    del domain_combination_sub[combination]
                    domain_combination_sub[f"{combination},{domain_id}"] = 1
                continue

            for combination in sorted(domain_combination_sub.keys()):
                del domain_combination_sub[combination]
                for domain in domains:
                    domain_id = split_domain_num(domain)
                    domain_combination_sub[f"{combination},{domain_id}"] = 1

        for combination in sorted(domain_combination_sub.keys()):
            domain_combination[combination] = 1

    return domain_combination


def load_rule_records(rule_file):
    rule_obj = {}

    with open(rule_file, 'r') as fh:
        rule_id = name = family = required = auxiliary = forbidden = rule_type = desc = ''

        for line in fh:
            line = line.strip()
            if line.startswith('#') or not line:
                continue

            if line.startswith('//'):
                if not all([rule_id, name, family, required, auxiliary, forbidden, rule_type, desc]):
                    raise ValueError(f"Error: Undefined rule member {rule_id}")
                rule_obj[rule_id] = RuleDefinition(
                    id=rule_id,
                    name=name,
                    family=family,
                    required=parse_domain_rule(required),
                    auxiliary=parse_domain_rule(auxiliary),
                    forbidden=parse_domain_rule(forbidden),
                    type=rule_type,
                    desc=desc,
                )
                rule_id = name = family = required = auxiliary = forbidden = rule_type = desc = ''
            elif line.startswith('ID:'):
                rule_id = line.replace('ID:', '').strip()
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
                rule_type = line.replace('Type:', '').strip()
            elif line.startswith('Desc:'):
                desc = line.replace('Desc:', '').strip()

    return rule_obj


def build_domain_hit_collection(hits, scores=None):
    if isinstance(hits, DomainHitCollection):
        return hits

    if isinstance(hits, str):
        hit_values = tuple(hits.split('\t')) if hits else ()
    else:
        hit_values = tuple(hits)

    if scores is None:
        score_values = ()
    elif isinstance(scores, str):
        score_values = tuple(map(float, scores.split('\t'))) if scores else ()
    else:
        score_values = tuple(float(score) for score in scores)

    if len(hit_values) != len(score_values):
        raise ValueError("[ERR]hit num do not match score num")

    counts = {}
    scores_by_domain = {}
    for hit_id, score_value in zip(hit_values, score_values):
        counts[hit_id] = counts.get(hit_id, 0) + 1
        scores_by_domain.setdefault(hit_id, []).append(score_value)

    return DomainHitCollection(
        hits=hit_values,
        scores=score_values,
        counts=counts,
        scores_by_domain=scores_by_domain,
    )


def compare_array(hit_collection, required_domains):
    collection = build_domain_hit_collection(hit_collection)
    required_counts = {}
    for domain_id in required_domains:
        required_counts[domain_id] = required_counts.get(domain_id, 0) + 1

    match = 0
    match_score = 0

    for domain_id in sorted(required_counts.keys()):
        if domain_id in collection.counts and collection.counts[domain_id] >= required_counts[domain_id]:
            match += 1
            sorted_scores = sorted(collection.scores_by_domain[domain_id], reverse=True)
            needed = required_counts[domain_id]
            for score_value in sorted_scores:
                match_score += score_value
                needed -= 1
                if needed == 0:
                    break

    match_status = 0
    if match > 0:
        match_status = 1
    if match == len(required_counts):
        match_status = 2

    return match_status, match_score


def compare_rule(hmm_hit, hmm_score, hmm_hit_s, hmm_score_s, rule_pack):
    rule_id = 'NA'

    hit_collection = build_domain_hit_collection(hmm_hit, hmm_score)
    hit_collection_s = build_domain_hit_collection(hmm_hit_s, hmm_score_s)

    total_domain = 0
    total_score = 0

    for rid in sorted(rule_pack.keys()):
        required_h = rule_attr(rule_pack[rid], 'required')
        forbidden_h = rule_attr(rule_pack[rid], 'forbidden')

        f_status = 0
        if forbidden_h != 'NA':
            for forbidden in sorted(forbidden_h.keys()):
                match_status, match_score = compare_array(hit_collection, forbidden.split(','))
                if match_status > 0:
                    f_status = 1
        if f_status == 1:
            continue

        for required in sorted(required_h.keys()):
            required_domains = required.split(',')
            match_status, match_score = compare_array(hit_collection, required_domains)

            match_status_b, match_score_b = (0, 0)
            if rid in ['T0008', 'T0009', 'T0011', 'T0023']:
                match_status_b, match_score_b = compare_array(hit_collection_s, required_domains)

            domain_num = len(required_domains)
            if match_status == 2 or match_status_b == 2:
                if rid == 'T9999':
                    if total_domain == 0 and total_score == 0:
                        total_domain = domain_num
                        total_score = match_score
                        rule_id = rid
                    elif domain_num > total_domain and match_score > total_score:
                        total_domain = domain_num
                        total_score = match_score
                        rule_id = rid
                else:
                    if total_domain == 0 and total_score == 0:
                        total_domain = domain_num
                        total_score = match_score
                        rule_id = rid
                    elif domain_num >= total_domain and match_score > total_score:
                        total_domain = domain_num
                        total_score = match_score
                        rule_id = rid

    return rule_id
