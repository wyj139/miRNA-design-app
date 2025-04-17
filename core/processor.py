from .utils import save_common_and_off_target, save_to_text, load_from_text
from .build_kmer_dict import build_kmer_dictionary
from .find_targets import find_common_kmers_and_off_targets
from .config import OUTPUT_RESULT, OUTPUT_TEXT_FILE

class KmerProcessor:
    @staticmethod
    def build_dictionary(file_path):
        """构建并保存k-mer字典"""
        kmer_indices, gene_kmers = build_kmer_dictionary(file_path)
        if kmer_indices and gene_kmers:
            save_to_text(kmer_indices, gene_kmers, OUTPUT_TEXT_FILE)
            return {
                'sequence_count': len(gene_kmers),
                'kmer_count': len(kmer_indices)
            }
        return None

    @staticmethod
    def load_dictionary():
        """加载已有字典"""
        return load_from_text(OUTPUT_TEXT_FILE)

    @staticmethod
    def process_gene_pair(gene1_id, gene2_id):
        """处理基因对，查找共同k-mer和脱靶基因"""
        kmer_indices, gene_kmers = load_from_text(OUTPUT_TEXT_FILE)
        if not (kmer_indices and gene_kmers):
            return None

        common_kmers, off_target_genes = find_common_kmers_and_off_targets(
            gene_kmers, gene1_id, gene2_id
        )

        if common_kmers is None:
            common_kmers = {}
        if off_target_genes is None:
            off_target_genes = {}
        
        

        save_common_and_off_target(common_kmers, off_target_genes, OUTPUT_RESULT)
        return common_kmers, off_target_genes