from Bio import SeqIO # type: ignore
from collections import defaultdict
from typing import List, Tuple, Dict

def reverse_complement(sequence):
    # 计算DNA序列的反向互补序列
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(sequence))

def read_fasta_file(file_path):
    # 读取FASTA文件并格式化基因ID，只保留相同gene_id中最长的对象
    try:
        sequences = {}
        for record in SeqIO.parse(file_path, "fasta"):
            parts = record.id.split('_')
            gene_id = '_'.join(parts[-2:])
            record.id = gene_id
            record.description = gene_id
            sequences[gene_id] = record
        return list(sequences.values())
    except FileNotFoundError:
        print(f"错误：找不到文件 {file_path}")
        return None
    except Exception as e:
        print(f"读取文件时发生错误：{str(e)}")
        return None 
def save_to_text(kmer_indices, gene_kmers, filename):
    """保存结果到文本文件"""
    try:
        if not kmer_indices or not gene_kmers:
            print(f"警告：空字典 - kmer_indices: {bool(kmer_indices)}, gene_kmers: {bool(gene_kmers)}")
            return False

        print(f"开始写入文件：{filename}")
        print(f"kmer_indices大小：{len(kmer_indices)}")
        print(f"gene_kmers大小：{len(gene_kmers)}")

        with open(filename, 'w', buffering=1) as f:  # 使用行缓冲
            # 保存kmer_indices
            print("写入kmer_indices...")
            f.write("7-mer to Gene Mapping:\n")
            for kmer, details in kmer_indices.items():
                f.write(f"{kmer}:\n")
                for detail in details:
                    f.write(f"    {detail}\n")
            f.write("\n")
            f.flush()
            print("kmer_indices写入完成")

            # 保存gene_kmers
            print("写入gene_kmers...")
            f.write("Gene to 7-mer Details:\n")
            for gene_id, kmers_info in gene_kmers.items():
                f.write(f"Gene ID: {gene_id}\n")
                for kmer, sequence, start_position in kmers_info:
                    f.write(f"  {kmer}: {sequence} ({start_position})\n")
                f.write("\n")
                f.flush()
            print("gene_kmers写入完成")

        # 验证文件内容
        print("验证文件内容...")
        with open(filename, 'r') as f:
            content = f.read()
            print(f"文件大小：{len(content)} 字节")
            print(f"文件是否包含Gene to 7-mer Details: {'Gene to 7-mer Details:' in content}")

        return True

    except Exception as e:
        print(f"保存文件时发生错误：{str(e)}")
        return False

def save_common_and_off_target(common_kmers, off_target_genes, filename):
    """保存共同7-mer和脱靶基因结果到文本文件"""
    with open(filename, "w") as f:
        # 保存共同的 7-mer 信息
        f.write("Common 7-mers:\n")
        for kmer, details in common_kmers.items():
            f.write(f"{kmer}:\n")
            for detail in details:
                f.write(f"    {detail}\n")
            # 显示完整的脱靶信息
            off_target_count = len(off_target_genes.get(kmer, {}))
            f.write(f"    脱靶基因数量: {off_target_count}\n")
            f.write("\n")

        # 保存完整的脱靶序列信息
        f.write("Off-Target Genes:\n")
        for kmer, gene_info in off_target_genes.items():
            f.write(f"{kmer}:\n")
            for gene_id, sequences in gene_info.items():
                for seq, pos in sequences:
                    # 保存完整的序列信息，便于后续脱靶检测
                    f.write(f"    {gene_id}: {seq} ({pos})\n")
            f.write("\n")

def load_from_text(filename):
    """从文本文件中读取数据，并返回 kmer_indices 和 gene_kmers 字典"""
    kmer_indices = defaultdict(list)
    gene_kmers = defaultdict(set)

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"错误：文件 {filename} 未找到！")
        return None, None

    if not lines:
        print(f"错误：文件 {filename} 为空！")
        return None, None

    mode = None
    current_kmer = None
    current_gene_id = None

    for line in lines:
        line = line.strip()
        if not line:
            continue  # 跳过空行

        if line == "7-mer to Gene Mapping:":
            mode = "kmer_to_gene"
            continue
        elif line == "Gene to 7-mer Details:":
            mode = "gene_to_kmer"
            continue

        if mode == "kmer_to_gene":
            if line.endswith(":"):
                current_kmer = line[:-1]
                kmer_indices[current_kmer] = []
            else:
                kmer_indices[current_kmer].append(line.strip())

        elif mode == "gene_to_kmer":
            if line.startswith("Gene ID:"):
                current_gene_id = line.split(":")[1].strip()
                gene_kmers[current_gene_id] = set()
            else:
                try:
                    kmer_info = line.split(":")
                    kmer = kmer_info[0].strip()
                    sequence, position = kmer_info[1].strip().split(" (")
                    position = int(position.rstrip(")"))
                    gene_kmers[current_gene_id].add((kmer, sequence, position))
                except Exception as e:
                    print(f"解析错误：{line}，错误信息：{e}")

    if not kmer_indices and not gene_kmers:
        print("错误：解析失败，返回空字典！")
        return None, None

    return kmer_indices, gene_kmers

def highlight_sequence(sequence: str, regions: List[Tuple[int, int]]) -> List[Dict[str, str]]:
    """
    高亮标记序列中的关键区域。
    :param sequence: 输入的序列
    :param regions: 需要高亮的区域列表，每个区域为 (start, end)，1-based 索引
    :return: 包含高亮区域的列表，每个区域包含起止位置和对应的序列
    """
    highlights = []
    for start, end in regions:
        if start <= len(sequence) and end <= len(sequence):
            highlights.append({
                "region": f"{start}-{end}",
                "sequence": sequence[start - 1:end]
            })
    return highlights