from Bio import SeqIO
from Bio.Seq import Seq


def is_nucleotide(sequence):
    nucleotide_set = set("ACGTUacgtu")
    return set(sequence).issubset(nucleotide_set)


def seq_to_hash(input_file):
    seq_info = {}
    with open(input_file) as input_handle:
        for inseq in SeqIO.parse(input_handle, 'fasta'):
            alphabet = 'protein'
            if is_nucleotide(str(inseq.seq)):
                alphabet = 'nucleotide'
            seq_info[inseq.id] = {
                'alphabet': alphabet,
                'seq': str(inseq.seq),
            }
    return seq_info


def translate_three_frames(dna_sequence):
    frames = []
    for i in range(3):
        seq_mod = dna_sequence[i:]
        if len(seq_mod) % 3 != 0:
            seq_mod += "N" * (3 - len(seq_mod) % 3)
        frames.append(seq_mod.translate())
    return frames


def six_frame_translation(dna_sequence):
    return translate_three_frames(dna_sequence) + translate_three_frames(dna_sequence.reverse_complement())


def translate_frames(dna_sequence, frame_mode):
    if frame_mode == '3F':
        return translate_three_frames(dna_sequence)
    if frame_mode == '3R':
        return translate_three_frames(dna_sequence.reverse_complement())
    return six_frame_translation(dna_sequence)


def prepare_input_sequences(input_file, input_protein_f, frame_mode):
    seq_info = seq_to_hash(input_file)
    input_count = len(seq_info)
    protein_num, nucleotide_num = 0, 0

    with open(input_protein_f, 'w') as outp:
        for seq_id in sorted(seq_info.keys()):
            if seq_info[seq_id]['alphabet'] == 'protein':
                protein_num += 1
                outp.write(f">{seq_id}\n{seq_info[seq_id]['seq']}\n")
                continue

            nucleotide_num += 1
            translated_frames = translate_frames(Seq(seq_info[seq_id]['seq']), frame_mode)
            for frame_index, translated_frame in enumerate(translated_frames, start=1):
                translated_id = f"{seq_id}-{frame_index}"
                outp.write(f">{translated_id}\n{translated_frame}\n")
                seq_info[translated_id] = {'alphabet': 'protein', 'seq': str(translated_frame)}

    return seq_info, input_count, protein_num, nucleotide_num
