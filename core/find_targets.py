from collections import defaultdict
from typing import Dict, List, Tuple, Set

def find_common_kmers_and_off_targets(gene_kmers, gene1_id, gene2_id):
    """
    查找两个基因之间的共同 7-mer 和脱靶基因，并返回更详细的脱靶信息
    :return: (common_kmers, off_targets)
        common_kmers: Dict[str, List[str]] - 共同 7-mer 的详细信息
        off_targets: Dict[str, Dict[str, List[Tuple[str, int]]]] - 每个 7-mer 对应的脱靶基因信息
    """
    if gene1_id not in gene_kmers or gene2_id not in gene_kmers:
        print("can't find")
        return None, None

    # 获取两个基因的 7-mer 集合
    gene1_kmers = {info[0] for info in gene_kmers[gene1_id]}
    gene2_kmers = {info[0] for info in gene_kmers[gene2_id]}

    # 找到共同的 7-mer
    common_kmers_set = gene1_kmers & gene2_kmers

    # 构建包含详细信息的字典
    common_kmers = {}
    off_targets = {}
    
    for kmer in common_kmers_set:
        # 处理共同 7-mer 信息
        kmer_details = []
        for gene_id in [gene1_id, gene2_id]:
            for kmer_info in gene_kmers[gene_id]:
                if kmer_info[0] == kmer:
                    formatted_entry = f"{gene_id}:{kmer_info[1]} ({kmer_info[2]})"
                    kmer_details.append(formatted_entry)
        common_kmers[kmer] = kmer_details

        # 处理脱靶基因信息
        off_targets[kmer] = defaultdict(list)
        for gene_id, kmers_info in gene_kmers.items():
            if gene_id not in (gene1_id, gene2_id):
                for kmer_info in kmers_info:
                    if kmer_info[0] == kmer:
                        # 存储更详细的脱靶信息：完整序列和位置
                        off_targets[kmer][gene_id].append((kmer_info[1], kmer_info[2]))

    return common_kmers, off_targets