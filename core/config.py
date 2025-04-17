# 文件路径配置
INPUT_FILE = "data/NCBI.RefSeq.Curated.3UTR.fasta"
OUTPUT_TEXT_FILE = "data/kmer_dictionary.txt"
OUTPUT_RESULT = "data/common_and_off_target_results.txt"

# 序列参数配置
KMER_LENGTH = 7
SEQUENCE_WINDOW = 30  # 提取的序列窗口长度

# 7-mer 规则配置
VALID_FIRST_POSITIONS = {'A', 'T'}  # 第1、3、6位允许的碱基
VALID_LAST_POSITION = {'G', 'C'}    # 第7位允许的碱基 
VALID_BASES = {'A', 'T', 'C', 'G'}  # 允许的所有碱基

# miRNA设计相关配置
MAX_MIRNA_LEN = 30  # 最大长度
GC_CONTENT_RANGE = (30.0, 70.0)  # GC含量范围（百分比）

# 脱靶判定配置
MIN_OFF_TARGET_MATCHES = 15  # 判定为脱靶所需的最小匹配碱基数

# 关键区域配置（1-based索引）
MIRNA_REGIONS = {        
    "supplementary": (9, 13),   # 10–14位
    "anchor": (16, 16),         # 17位
    "tail": (18, 20)            # 19–21位
}

# 高亮显示配置
HIGHLIGHT_POSITIONS = MIRNA_REGIONS  # 使用相同的区域定义
