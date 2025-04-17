from typing import Dict, List, Tuple, Set, Any
from .config import HIGHLIGHT_POSITIONS, MAX_MIRNA_LEN, VALID_BASES, MIN_OFF_TARGET_MATCHES

class InteractiveDesigner:
    def __init__(self, selected_kmer: str, common_kmers: Dict, off_targets: Dict):
        """
        初始化交互设计器
        :param selected_kmer: 用户选定的 7-mer 序列
        :param common_kmers: 共同 7-mer 信息
        :param off_targets: 脱靶基因信息
        """
        self.selected_kmer = selected_kmer
        self.current_sequence = ""  # 用户添加的序列
        self.common_kmers = common_kmers
        self.off_targets = off_targets
        
        # 存储该7mer对应的所有脱靶序列信息
        self.off_target_sequences = self.get_off_target_sequences(selected_kmer)
        
        self.sequence_changed = False

    def get_off_target_sequences(self, kmer: str) -> Dict[str, List[Tuple[str, int]]]:
        """
        获取指定 7-mer 的所有脱靶序列信息
        :return: {gene_id: [(sequence, position), ...]}
        """
        return self.off_targets.get(kmer, {})

    def add_base(self, base: str) -> bool:
        """
        添加新碱基并更新脱靶信息
        返回是否添加成功
        """
        if not self._validate_base(base):
            return False
            
        if len(self.get_full_sequence()) >= MAX_MIRNA_LEN:
            return False
            
        self.current_sequence += base
        self.sequence_changed = True
        return True

    def _validate_base(self, base: str) -> bool:
        """验证碱基是否合法"""
        return base in VALID_BASES

    def get_full_sequence(self) -> str:
        """
        获取完整序列（7-mer + 当前序列）
        返回：先是7-mer，后面是用户添加的碱基
        """
        return self.selected_kmer + self.current_sequence

    def get_gc_content(self) -> float:
        """计算完整序列的GC含量"""
        full_sequence = self.get_full_sequence()
        if not full_sequence:
            return 0.0
        gc_count = sum(1 for base in full_sequence if base in "GC")
        return (gc_count / len(full_sequence)) * 100

    def get_off_targets(self) -> Tuple[int, List[str]]:
        """
        获取当前脱靶基因信息
        通过比较当前完整序列与脱靶序列的匹配程度来确定脱靶基因
        需同时满足：
        1. 至少有两组连续匹配（允许±1位置偏移）
        2. 总匹配碱基数不少于15个
        """
        current_off_targets = set()
        full_sequence = self.get_full_sequence()
        
        # 遍历所有脱靶基因及其序列
        for gene_id, sequences in self.off_target_sequences.items():
            for seq_info in sequences:
                target_seq = seq_info[0]
                
                # 首先检查总体匹配碱基数
                total_matches = sum(1 for i in range(min(len(full_sequence), len(target_seq)))
                                  if full_sequence[i] == target_seq[i])
                
                # 如果总匹配碱基数达到要求，再检查连续匹配组
                if total_matches >= MIN_OFF_TARGET_MATCHES:
                    # 只比较7-mer后的序列部分
                    added_sequence = full_sequence[len(self.selected_kmer):]
                    target_added = target_seq[len(self.selected_kmer):]
                    
                    if added_sequence:  # 确保有新增序列
                        # 查找连续匹配组
                        match_groups = self._find_match_groups_with_shift(added_sequence, target_added)
                        
                        # 如果同时满足两个条件：总匹配数>=15 且 有至少两组连续匹配
                        if len(match_groups) >= 2:
                            current_off_targets.add(gene_id)
                            break

        return len(current_off_targets), list(sorted(current_off_targets))

    def _find_match_groups_with_shift(self, seq1: str, seq2: str) -> List[Tuple[int, int, int]]:
        """
        查找两个序列中所有连续三个碱基匹配的组，允许±1的位置偏移
        :return: [(seq1_start, seq2_start, length), ...]
        """
        match_groups = []
        seq1_len = len(seq1)
        seq2_len = len(seq2)
        
        for seq1_pos in range(seq1_len - 2):  # 至少需要3个碱基
            # 在目标序列中检查前后1个位置
            for shift in [-1, 0, 1]:
                seq2_pos = seq1_pos + shift
                if seq2_pos < 0 or seq2_pos + 2 >= seq2_len:
                    continue
                
                # 检查连续3个碱基是否匹配
                if all(seq1[seq1_pos + i] == seq2[seq2_pos + i] for i in range(3)):
                    match_groups.append((seq1_pos, seq2_pos, 3))

        # 合并相邻或重叠的匹配组
        return self._merge_adjacent_groups(match_groups)

    def _merge_adjacent_groups(self, groups: List[Tuple[int, int, int]]) -> List[Tuple[int, int, int]]:
        """
        合并相邻或重叠的匹配组
        """
        if not groups:
            return []
            
        # 按seq1的起始位置排序
        groups.sort()
        merged = []
        current = list(groups[0])  # [seq1_start, seq2_start, length]
        
        for next_group in groups[1:]:
            # 如果当前组和下一组的位置差不超过1
            if (next_group[0] - (current[0] + current[2])) <= 1 and \
               (next_group[1] - (current[1] + current[2])) <= 1:
                # 扩展当前匹配组
                current[2] = next_group[0] + 3 - current[0]
            else:
                # 保存当前组并开始新的组
                merged.append(tuple(current))
                current = list(next_group)
                
        merged.append(tuple(current))
        return merged

    def get_highlighted_sequence(self) -> List[Dict[str, Any]]:
        """返回格式化的高亮信息，标记关键区域"""
        highlights = []
        full_sequence = self.get_full_sequence()  # 获取完整序列
        for region, (start, end) in HIGHLIGHT_POSITIONS.items():
            if end <= len(full_sequence):  # 使用完整序列长度
                highlights.append({
                    "region": (start, end),
                    "sequence": full_sequence[start - 1:end]  # 使用完整序列
                })
        return highlights

    def reset_sequence(self):
        """重置序列和脱靶基因集合"""
        self.current_sequence = ""
        self.off_target_sequences = self.get_off_target_sequences(self.selected_kmer)
        self.sequence_changed = True

    def get_off_target_details(self, gene_ids: List[str]) -> Dict[str, Dict]:
        """
        获取脱靶基因的详细匹配信息
        :param gene_ids: 脱靶基因ID列表
        :return: {
            gene_id: {
                "sequence": 完整匹配序列,
                "matches": 匹配位置信息
            }
        }
        """
        details = {}
        full_sequence = self.get_full_sequence()
        added_sequence = full_sequence[len(self.selected_kmer):]

        for gene_id in gene_ids:
            if gene_id in self.off_target_sequences:
                for seq_info in self.off_target_sequences[gene_id]:
                    target_seq = seq_info[0]
                    target_added = target_seq[len(self.selected_kmer):]
                    
                    # 查找匹配组
                    match_groups = self._find_match_groups_with_shift(
                        added_sequence, target_added
                    )
                    
                    if match_groups:
                        details[gene_id] = {
                            "sequence": target_seq,
                            "matches": [
                                f"{g[0]+len(self.selected_kmer)}-{g[0]+g[2]+len(self.selected_kmer)}"
                                for g in match_groups
                            ]
                        }
                        break

        return details