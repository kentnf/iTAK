import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq

from itak.classify import aln_to_hash, identify_protein_kinases, itak_pk_classify
from itak.hmmscan import HmmscanAlignment, HmmscanSummaryHit
from itak.rules import build_domain_hit_collection, compare_array, compare_rule, load_rule_records
from itak.runtime import resolve_output_dir
from itak.sequences import six_frame_translation, translate_frames, translate_three_frames


class ITAKHelperTests(unittest.TestCase):
    def test_translate_frames_dispatches_to_expected_mode(self):
        sequence = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")

        self.assertEqual(
            [str(frame) for frame in translate_frames(sequence, "6")],
            [str(frame) for frame in six_frame_translation(sequence)],
        )
        self.assertEqual(
            [str(frame) for frame in translate_frames(sequence, "3F")],
            [str(frame) for frame in translate_three_frames(sequence)],
        )
        self.assertEqual(
            [str(frame) for frame in translate_frames(sequence, "3R")],
            [str(frame) for frame in translate_three_frames(sequence.reverse_complement())],
        )

    def test_resolve_output_dir_keeps_multi_input_outputs_separate(self):
        self.assertEqual(
            resolve_output_dir("input_a.fa", "/tmp/out", True),
            "/tmp/out/input_a.fa_output",
        )
        self.assertEqual(
            resolve_output_dir("input_a.fa", "/tmp/out", False),
            "/tmp/out",
        )

    def test_identify_protein_kinases_applies_ga_cutoffs(self):
        hmmscan_hit = [
            HmmscanSummaryHit("gene1", "PF00069.29", "120.0", "1e-40"),
            HmmscanSummaryHit("gene2", "PF07714.21", "95.0", "1e-25"),
            HmmscanSummaryHit("gene3", "PF00069.29", "10.0", "1e-3"),
            HmmscanSummaryHit("gene4", "PF99999.1", "200.0", "1e-60"),
        ]
        ga_cutoff = {"PF00069": 50.0, "PF07714": 90.0}

        self.assertEqual(
            identify_protein_kinases(hmmscan_hit, ga_cutoff),
            {"gene1": 1, "gene2": 1},
        )

    def test_identify_protein_kinases_requires_ga_for_kinase_domain(self):
        with self.assertRaises(Exception):
            identify_protein_kinases(
                [HmmscanSummaryHit("gene1", "PF00069.29", "120.0", "1e-40")],
                {},
            )

    def test_aln_to_hash_accepts_structured_alignments(self):
        alignments = [
            HmmscanAlignment(
                query_name="gene1",
                hit_id="PF00069.29",
                query_start="19",
                query_end="274",
                hit_start="1",
                hit_end="263",
                query_seq="274",
                alignment="",
                hit_seq="263",
                score="247.2",
                evalue="1.2e-75",
                description="Protein kinase domain",
                query_length="448",
            )
        ]

        self.assertEqual(
            aln_to_hash(alignments, {"PF00069": 100.0}),
            {
                "gene1": "gene1\tPF00069.29\t19\t274\t1\t263\t274\t\t263\t247.2\t1.2e-75\tProtein kinase domain\t448\n"
            },
        )

    def test_compare_array_uses_highest_scores_for_required_duplicates(self):
        collection = build_domain_hit_collection(
            ["PF1", "PF1", "PF1", "PF2"],
            [10.0, 50.0, 20.0, 5.0],
        )

        match_status, match_score = compare_array(collection, ["PF1", "PF1"])

        self.assertEqual(match_status, 2)
        self.assertEqual(match_score, 70.0)

    def test_compare_rule_respects_forbidden_domains_and_prefers_higher_score(self):
        rule_pack = {
            "T0001": {
                "name": "Rule1",
                "family": "Fam1",
                "required": {"PF1,PF2": 1},
                "auxiliary": {},
                "forbidden": {},
                "type": "TF",
                "desc": "rule1",
            },
            "T0002": {
                "name": "Rule2",
                "family": "Fam2",
                "required": {"PF1,PF2": 1},
                "auxiliary": {},
                "forbidden": {"PFX": 1},
                "type": "TF",
                "desc": "rule2",
            },
            "T9999": {
                "name": "Others",
                "family": "Others",
                "required": {"PF1": 1},
                "auxiliary": {},
                "forbidden": {},
                "type": "TR",
                "desc": "fallback",
            },
        }

        hits = build_domain_hit_collection(["PF1", "PF2", "PFX"], [10.0, 20.0, 30.0])
        best_hits = build_domain_hit_collection(["PF1", "PF2"], [10.0, 20.0])
        self.assertEqual(compare_rule(hits, None, best_hits, None, rule_pack), "T0001")

        stronger_hits = build_domain_hit_collection(["PF1"], [99.0])
        self.assertEqual(compare_rule(stronger_hits, None, stronger_hits, None, {"T9999": rule_pack["T9999"]}), "T9999")

    def test_itak_pk_classify_accepts_structured_alignments(self):
        alignments = [
            HmmscanAlignment(
                query_name="gene1",
                hit_id="PK_A",
                query_start="1",
                query_end="10",
                hit_start="1",
                hit_end="10",
                query_seq="10",
                alignment="",
                hit_seq="10",
                score="100.0",
                evalue="1e-20",
                description="desc",
                query_length="100",
            ),
            HmmscanAlignment(
                query_name="gene1",
                hit_id="PK_B",
                query_start="1",
                query_end="10",
                hit_start="1",
                hit_end="10",
                query_seq="10",
                alignment="",
                hit_seq="10",
                score="120.0",
                evalue="1e-25",
                description="desc",
                query_length="100",
            ),
        ]

        hit, best_align = itak_pk_classify(alignments, {"gene1": 1, "gene2": 1}, "Group-other")

        self.assertEqual(hit, {"gene1": "PK_B", "gene2": "Group-other"})
        self.assertEqual(
            best_align["gene1"],
            "gene1\tPK_B\t1\t10\t1\t10\t10\t\t10\t120.0\t1e-25\tdesc\t100\n",
        )
        self.assertEqual(best_align["gene2"], "")

    def test_load_rule_records_can_serialize_legacy_dict(self):
        rule_text = """\
ID:T0001
Name:TestRule
Family:TestFam
Required:PF1#1:PF2#1
Auxiliary:NA
Forbidden:PF3#1
Type:TF
Desc:demo
//
"""

        with tempfile.TemporaryDirectory() as temp_dir:
            rule_path = Path(temp_dir) / "rule.txt"
            rule_path.write_text(rule_text)

            records = load_rule_records(rule_path)

        self.assertEqual(set(records.keys()), {"T0001"})
        self.assertEqual(records["T0001"].name, "TestRule")
        self.assertEqual(records["T0001"].family, "TestFam")
        self.assertEqual(records["T0001"].required, {"PF1": 1, "PF2": 1})
        self.assertEqual(records["T0001"].forbidden, {"PF3": 1})
        self.assertEqual(
            records["T0001"].to_legacy_dict(),
            {
                "name": "TestRule",
                "family": "TestFam",
                "required": {"PF1": 1, "PF2": 1},
                "auxiliary": {},
                "forbidden": {"PF3": 1},
                "type": "TF",
                "desc": "demo",
            },
        )


if __name__ == "__main__":
    unittest.main()
