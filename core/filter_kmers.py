from .config import VALID_FIRST_POSITIONS, VALID_LAST_POSITION

def is_valid_kmer(kmer):
    """
    检查 7-mer 是否符合规则：
    - 第 1、3、6 位必须是 A 或 T
    - 第 7 位必须是 G 或 C
    """
    return (kmer[0] in VALID_FIRST_POSITIONS and 
            kmer[2] in VALID_FIRST_POSITIONS and 
            kmer[5] in VALID_FIRST_POSITIONS and 
            kmer[6] in VALID_LAST_POSITION)

def filter_kmers(kmer_dict):
    """过滤符合规则的 7-mer"""
    return {kmer: details for kmer, details in kmer_dict.items() 
            if is_valid_kmer(kmer)} 