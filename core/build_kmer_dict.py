from collections import defaultdict
from .utils import reverse_complement, read_fasta_file
from .filter_kmers import is_valid_kmer
from .config import KMER_LENGTH, SEQUENCE_WINDOW


def build_kmer_dictionary(file_path, progress_callback=None):
    """构建 7-mer 字典"""
    
    sequences = read_fasta_file(file_path)
    if not sequences:
        print("fail")
        return None, None
    print(f"共处理{len(sequences)}条序列")
    
    kmer_indices = defaultdict(list)
    gene_kmers = defaultdict(set)
    total = len(sequences)
    
    for i, seq_record in enumerate(sequences):
        sequence = str(seq_record.seq)
        seq_id = seq_record.id
        sequence_length = len(sequence)
        
        rev_comp_sequence = reverse_complement(sequence)

        for i in range(len(rev_comp_sequence) - SEQUENCE_WINDOW + 1):
            kmer = rev_comp_sequence[i:i + KMER_LENGTH]
            if is_valid_kmer(kmer):
                original_position = sequence_length - (i + KMER_LENGTH)
                kmer_indices[kmer].append(
                    f"{seq_id}:{rev_comp_sequence[i:i + SEQUENCE_WINDOW]} ({original_position})"
                )
                gene_kmers[seq_id].add(
                    (kmer, rev_comp_sequence[i:i + SEQUENCE_WINDOW], original_position)
                )
        
        if progress_callback:
            progress = (i + 1) / total * 100
            progress_callback(progress)

    return kmer_indices, gene_kmers
