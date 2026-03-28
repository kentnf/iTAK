import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

import itak.sequences as itak_sequences


class ITAKSequenceTests(unittest.TestCase):
    def test_prepare_input_sequences_translates_nucleotide_inputs_in_six_frames(self):
        fasta_text = ">protein1\nMSTNPKPQR\n>mrna1\nATGGCCATT\n"

        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "input.fa"
            output_path = Path(temp_dir) / "protein.fa"
            input_path.write_text(fasta_text)

            seq_info, input_count, protein_num, nucleotide_num = itak_sequences.prepare_input_sequences(
                input_path,
                output_path,
                "6",
            )

            with open(output_path) as output_handle:
                written_ids = [record.id for record in SeqIO.parse(output_handle, "fasta")]

        expected_frames = [
            str(frame)
            for frame in itak_sequences.six_frame_translation(Seq("ATGGCCATT"))
        ]

        self.assertEqual(input_count, 2)
        self.assertEqual(protein_num, 1)
        self.assertEqual(nucleotide_num, 1)
        self.assertEqual(written_ids, ["mrna1-1", "mrna1-2", "mrna1-3", "mrna1-4", "mrna1-5", "mrna1-6", "protein1"])
        self.assertEqual(
            [seq_info[f"mrna1-{index}"]["seq"] for index in range(1, 7)],
            expected_frames,
        )

    def test_prepare_input_sequences_respects_reverse_three_frame_mode(self):
        fasta_text = ">mrna1\nATGGCCATT\n"

        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "input.fa"
            output_path = Path(temp_dir) / "protein.fa"
            input_path.write_text(fasta_text)

            seq_info, input_count, protein_num, nucleotide_num = itak_sequences.prepare_input_sequences(
                input_path,
                output_path,
                "3R",
            )

            with open(output_path) as output_handle:
                written_ids = [record.id for record in SeqIO.parse(output_handle, "fasta")]

        expected_frames = [
            str(frame)
            for frame in itak_sequences.translate_three_frames(Seq("ATGGCCATT").reverse_complement())
        ]

        self.assertEqual(input_count, 1)
        self.assertEqual(protein_num, 0)
        self.assertEqual(nucleotide_num, 1)
        self.assertEqual(written_ids, ["mrna1-1", "mrna1-2", "mrna1-3"])
        self.assertEqual(
            [seq_info[f"mrna1-{index}"]["seq"] for index in range(1, 4)],
            expected_frames,
        )


if __name__ == "__main__":
    unittest.main()
