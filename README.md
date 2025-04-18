# miRNA-design-app

## 项目简介

`miRNA-design-app` 是一个用于设计和分析 miRNA 的工具，支持从基因序列中提取 k-mer（7-mer），筛选共同序列，检测脱靶基因，并提供交互式 miRNA 设计功能。该工具结合了生物信息学的基本算法和用户友好的图形界面，适用于 miRNA 相关的研究和设计任务。

---

## 功能特性

1. **k-mer 字典构建**  
   从输入的基因序列文件（FASTA 格式）中提取 7-mer，并构建基因与 k-mer 的映射关系。

2. **共同 k-mer 和脱靶基因分析**  
   支持基因对的共同 k-mer 筛选，并检测脱靶基因及其详细信息。

3. **交互式 miRNA 设计**  
   提供交互界面，用户可以基于选定的 7-mer 添加碱基，实时查看 GC 含量、脱靶基因数量及匹配信息。

4. **结果保存与加载**  
   支持将分析结果保存到文本文件，并可加载已有的字典文件进行后续操作。

5. **图形用户界面 (GUI)**  
   提供基于 Tkinter 的用户界面，方便用户操作和查看结果。

---

## 文件结构
miRNA-design-app/ ├── core/ # 核心功能模块 │ ├── __init__.py # 模块初始化 │ ├── build_kmer_dict.py # k-mer 字典构建 │ ├── config.py # 配置文件 │ ├── filter_kmers.py # k-mer 过滤规则 │ ├── find_targets.py # 共同 k-mer 和脱靶基因分析 │ ├── interactive_designer.py # 交互式设计器 │ ├── processor.py # 数据处理器 │ └── utils.py # 工具函数 ├── data/ # 数据文件 │ ├── NCBI.RefSeq.Curated.3UTR.fasta # 示例输入文件 │ ├── kmer_dictionary.txt # k-mer 字典文件 │ ├── common_and_off_target_results.txt # 分析结果文件 ├── GUI/ # 图形用户界面 │ ├── app.py # 主界面逻辑 ├── tests/ # 单元测试 │ ├── test_build_kmer_dict.py # 测试 k-mer 字典构建 │ ├── test_find_targets.py # 测试共同 k-mer 和脱靶基因分析 │ ├── test_interactive_designer.py # 测试交互式设计器 │ ├── test_processor.py # 测试数据处理器 │ ├── test_utils.py # 测试工具函数 ├── main.py # 程序入口 ├── pytest.ini # pytest 配置 ├── README.md # 项目说明文档 └── .gitignore # Git 忽略文件

---

## 安装与运行

### 环境要求

- Python 3.8 或更高版本
- 必要的依赖库：
  - `biopython`
  - `tkinter`

### 安装依赖

使用以下命令安装所需依赖：

```bash
pip install biopython
```

---

## 运行程序
启动图形界面：

按照界面提示选择操作方式：

使用现有字典
重新构建字典
筛选共同序列
交互式设计 miRNA

---

## 使用说明
1. 构建 k-mer 字典
选择一个 FASTA 文件作为输入。
点击“创建字典”按钮，程序将提取 7-mer 并生成字典文件。
2. 筛选共同序列
输入两个基因 ID。
点击“筛选序列”按钮，程序将分析共同的 7-mer 和脱靶基因信息。
3. 交互式 miRNA 设计
选择一个 7-mer。
添加碱基并实时查看 GC 含量和脱靶基因信息。
支持高亮显示关键区域。

---

## 配置说明
配置文件位于 core/config.py，主要参数包括：

文件路径配置

INPUT_FILE：输入文件路径
OUTPUT_TEXT_FILE：k-mer 字典输出路径
OUTPUT_RESULT：分析结果输出路径
序列参数

KMER_LENGTH：k-mer 长度（默认为 7）
SEQUENCE_WINDOW：序列窗口长度
miRNA 设计参数

MAX_MIRNA_LEN：miRNA 最大长度
GC_CONTENT_RANGE：GC 含量范围
脱靶判定参数

MIN_OFF_TARGET_MATCHES：脱靶判定所需的最小匹配碱基数

---

## 测试
项目包含单元测试，位于 tests/ 目录。运行以下命令执行测试：

---

## 示例数据
示例数据位于 data/ 目录，包括：

NCBI.RefSeq.Curated.3UTR.fasta：FASTA 格式的基因序列文件
kmer_dictionary.txt：生成的 k-mer 字典
common_and_off_target_results.txt：共同 k-mer 和脱靶基因分析结果

---

## 贡献
欢迎提交 Issue 或 Pull Request 来改进本项目。

---
