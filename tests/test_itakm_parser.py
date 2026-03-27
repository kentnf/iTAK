import tempfile
import textwrap
import unittest
from pathlib import Path

import itakm


class ITAKMParserTests(unittest.TestCase):
    def test_parse_align_records_returns_structured_alignment(self):
        hsp_info = textwrap.dedent(
            """\
            >> PF00069.29 Protein kinase domain
               1 ! 247.2 0.0 1.2e-75 1.2e-75 1 263 .. 19 274 ..
            == domain 1
              PF00069.29 1 YELLEKLGSGSFGKV 16
                           ye++++lG+Gsf+kV
              gene1 19 YEMGRTLGEGSFAKV 34
            """
        )

        alignments = itakm.parse_align_records(hsp_info, "gene1", "448")

        self.assertEqual(len(alignments), 1)
        alignment = alignments[0]
        self.assertEqual(alignment.query_name, "gene1")
        self.assertEqual(alignment.hit_id, "PF00069.29")
        self.assertEqual(alignment.query_start, "19")
        self.assertEqual(alignment.query_end, "274")
        self.assertEqual(alignment.hit_start, "1")
        self.assertEqual(alignment.hit_end, "263")
        self.assertEqual(alignment.query_seq, "34")
        self.assertEqual(alignment.alignment, "")
        self.assertEqual(alignment.hit_seq, "16")
        self.assertEqual(alignment.score, "247.2")
        self.assertEqual(alignment.evalue, "1.2e-75")
        self.assertEqual(alignment.description, "Protein kinase domain")

    def test_parse_hmmscan_records_can_serialize_legacy_output(self):
        hmmscan_text = textwrap.dedent(
            """\
            Query: gene1 [L=448]
            >> PF00069.29 Protein kinase domain
               1 ! 247.2 0.0 1.2e-75 1.2e-75 1 263 .. 19 274 ..
            == domain 1
              PF00069.29 1 YELLEKLGSGSFGKV 16
                           ye++++lG+Gsf+kV
              gene1 19 YEMGRTLGEGSFAKV 34
            //
            """
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            result_path = Path(temp_dir) / "hmmscan.txt"
            result_path.write_text(hmmscan_text)

            parse_result = itakm.parse_hmmscan_records(result_path)
            legacy_summary, legacy_alignment, legacy_best = parse_result.to_legacy_tuple()

        self.assertEqual(len(parse_result.summary_hits), 1)
        self.assertEqual(len(parse_result.alignments), 1)
        self.assertEqual(len(parse_result.best_domain_hits), 0)
        self.assertEqual(
            legacy_summary,
            "gene1\tPF00069.29\t247.2\t1.2e-75\n",
        )
        self.assertEqual(
            legacy_alignment,
            "gene1\tPF00069.29\t19\t274\t1\t263\t34\t\t16\t247.2\t1.2e-75\tProtein kinase domain\t448\n",
        )
        self.assertEqual(legacy_best, "")

    def test_parse_hmmscan_records_handles_no_hits(self):
        hmmscan_text = textwrap.dedent(
            """\
            Query: gene2 [L=120]
            No hits detected that satisfy reporting thresholds
            //
            """
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            result_path = Path(temp_dir) / "hmmscan.txt"
            result_path.write_text(hmmscan_text)

            parse_result = itakm.parse_hmmscan_records(result_path)

        self.assertEqual(parse_result.summary_hits, [])
        self.assertEqual(parse_result.alignments, [])
        self.assertEqual(parse_result.best_domain_hits, [])


if __name__ == "__main__":
    unittest.main()
