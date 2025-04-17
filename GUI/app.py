import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from core.processor import KmerProcessor
from core.config import OUTPUT_TEXT_FILE
from core.interactive_designer import InteractiveDesigner
from core.config import HIGHLIGHT_POSITIONS

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("miRNA k-mer 处理")
        self.processor = KmerProcessor()
        self.designer = None  # InteractiveDesigner 实例
        
        # 初始化存储共同 7-mer 和脱靶信息的属性
        self.common_kmers = {}
        self.off_targets = {}
        
        self.setup_gui()

    def setup_gui(self):
        """设置GUI组件"""
        self.frame = tk.Frame(self.root, padx=20, pady=20)
        self.frame.pack(expand=True, fill='both')

        # 选项标签
        tk.Label(self.frame, text="请选择操作方式：", font=('Arial', 12)).pack(pady=10)

        # 主要按钮
        self.create_main_buttons()

        # 初始化其他组件
        self.setup_file_selection()
        self.setup_result_display()
        self.create_gene_filter_frame()

        # 添加交互设计模块
        self.setup_interactive_designer()

        self.file_path = ""

    def create_main_buttons(self):
        """创建主要按钮"""
        for text, command in [
            ("使用现有字典", self.use_existing_dict),
            ("重新构建字典", self.show_file_selection)
        ]:
            tk.Button(
                self.frame,
                text=text,
                command=command,
                width=20,
                height=2
            ).pack(pady=10)

    def setup_file_selection(self):
        """设置文件选择相关组件"""
        self.entry_var = tk.StringVar()
        self.entry = tk.Entry(self.frame, textvariable=self.entry_var, width=50)
        self.file_select_button = tk.Button(
            self.frame, 
            text="选择FASTA文件", 
            command=self.select_file
        )
        self.build_dict_button = tk.Button(
            self.frame, 
            text="创建字典", 
            command=self.build_dict
        )

    def setup_result_display(self):
        """设置结果显示区域"""
        self.result_text = tk.Text(self.frame, height=20, width=80)
        self.result_text.pack(pady=10)
        self.result_text.pack_forget()

    def create_gene_filter_frame(self):
        """创建基因筛选框架"""
        self.filter_frame = tk.Frame(self.frame)
        self.filter_frame.pack(pady=10)
        
        # 基因输入框
        for i in range(1, 3):
            tk.Label(self.filter_frame, text=f"输入基因{i}:").pack(side=tk.LEFT, padx=5)
            entry = tk.Entry(self.filter_frame, width=15)
            entry.pack(side=tk.LEFT, padx=5)
            setattr(self, f'gene{i}_entry', entry)

        # 筛选按钮
        tk.Button(
            self.filter_frame,
            text="筛选序列",
            command=self.filter_sequences
        ).pack(side=tk.LEFT, padx=5)

    def setup_interactive_designer(self):
        """设置交互设计模块"""
        self.interactive_frame = tk.Frame(self.frame, padx=10, pady=10, relief=tk.GROOVE, borderwidth=2)
        self.interactive_frame.pack(pady=10, fill=tk.X)

        tk.Label(self.interactive_frame, text="交互式 miRNA 设计", font=('Arial', 12, 'bold')).pack(pady=5)

        # 创建一个包含输入框和清空按钮的框架
        kmer_frame = tk.Frame(self.interactive_frame)
        kmer_frame.pack(pady=5)
        
        tk.Label(kmer_frame, text="输入 7-mer：", font=('Arial', 10)).pack(side=tk.LEFT, padx=5)
        self.selected_kmer_entry = tk.Entry(kmer_frame, width=20)
        self.selected_kmer_entry.pack(side=tk.LEFT, padx=5)
        
        # 添加清空按钮
        tk.Button(
            kmer_frame,
            text="清空",
            command=self.clear_kmer_entry,
            width=5
        ).pack(side=tk.LEFT, padx=5)

        # 输入碱基按钮
        self.base_input_frame = tk.Frame(self.interactive_frame)
        self.base_input_frame.pack(pady=5)
        for base in ['A', 'T', 'C', 'G']:
            tk.Button(
                self.base_input_frame,
                text=base,
                command=lambda b=base: self.add_base_to_sequence(b),
                width=5
            ).pack(side=tk.LEFT, padx=5)

        # 显示区域
        self.sequence_label = tk.Label(self.interactive_frame, text="当前序列：", font=('Arial', 10))
        self.sequence_label.pack(pady=5)

        self.gc_content_label = tk.Label(self.interactive_frame, text="GC 含量：0.00%", font=('Arial', 10))
        self.gc_content_label.pack(pady=5)

        self.off_target_label = tk.Label(self.interactive_frame, text="脱靶基因数量：0", font=('Arial', 10))
        self.off_target_label.pack(pady=5)

        self.off_target_list = tk.Text(self.interactive_frame, height=5, width=50)
        self.off_target_list.pack(pady=5)

        # 重置按钮
        tk.Button(
            self.interactive_frame,
            text="重置序列",
            command=self.reset_sequence,
            width=15
        ).pack(pady=5)

    def clear_kmer_entry(self):
        """清空7-mer输入框并重置设计器"""
        self.selected_kmer_entry.delete(0, tk.END)
        self.designer = None
        self.update_interactive_display()
        # 重置显示
        self.sequence_label.config(text="当前序列：")
        self.gc_content_label.config(text="GC 含量：0.00%")
        self.off_target_label.config(text="脱靶基因数量：0")
        self.off_target_list.delete(1.0, tk.END)

    def add_base_to_sequence(self, base):
        """向当前序列添加碱基"""
        if not self.designer:
            # 确保用户选择了 7-mer
            selected_kmer = self.get_selected_kmer()
            if not selected_kmer:
                messagebox.showerror("错误", "请先选择一个 7-mer")
                return
                
            # 确保已经有共同 7-mer 和脱靶信息
            if not self.common_kmers or not self.off_targets:
                messagebox.showerror("错误", "请先进行序列筛选")
                return
                
            if selected_kmer not in self.common_kmers:
                messagebox.showerror("错误", "所选 7-mer 不在共同序列中")
                return

            # 初始化 InteractiveDesigner
            self.designer = InteractiveDesigner(
                selected_kmer,
                self.common_kmers,
                self.off_targets
            )

        # 添加碱基到当前序列
        self.designer.add_base(base)
        
        # 更新显示
        self.update_interactive_display()

    def get_selected_kmer(self):
        """获取用户选定的 7-mer"""
        # 假设在 GUI 中有一个输入框或下拉菜单供用户选择 7-mer
        # 这里以一个简单的输入框为例
        selected_kmer = getattr(self, "selected_kmer_entry", None)
        if selected_kmer:
            return selected_kmer.get().strip()
        return None

    def update_interactive_display(self):
        """更新交互设计模块的显示"""
        if not self.designer:
            return
        
        self._update_sequence_display()
        self._update_gc_content()
        self._update_off_target_info()

    def _update_sequence_display(self):
        """更新序列显示"""
        sequence = self.designer.get_full_sequence()
        highlighted_sequence = self._get_highlighted_sequence(sequence)
        kmer_len = len(self.designer.selected_kmer)
        
        display_text = (f"当前序列："
                       f"{highlighted_sequence[:kmer_len]}|"
                       f"{highlighted_sequence[kmer_len:]}")
        self.sequence_label.config(text=display_text)

    def _get_highlighted_sequence(self, sequence: str) -> str:
        """生成高亮显示的序列"""
        result = ""
        for i, base in enumerate(sequence):
            if self._is_position_highlighted(i + 1):
                result += f"[{base}]"
            else:
                result += base
        return result

    def _is_position_highlighted(self, position: int) -> bool:
        """检查位置是否需要高亮"""
        return any(start <= position <= end 
                  for start, end in HIGHLIGHT_POSITIONS.values())

    def _update_gc_content(self):
        """更新 GC 含量显示"""
        if self.designer:
            gc_content = self.designer.get_gc_content()
            self.gc_content_label.config(text=f"GC 含量：{gc_content:.2f}%")

    def _update_off_target_info(self):
        """更新脱靶基因信息显示"""
        if self.designer:
            count, genes = self.designer.get_off_targets()
            
            # 更新脱靶基因数量
            self.off_target_label.config(text=f"脱靶基因数量：{count}")
            
            # 更新脱靶基因列表，同时显示匹配序列信息
            self.off_target_list.delete(1.0, tk.END)
            if genes:
                # 获取详细的脱靶信息
                off_target_details = self.designer.get_off_target_details(genes)
                for gene_id, details in off_target_details.items():
                    match_info = details["matches"]
                    sequence = details["sequence"]
                    self.off_target_list.insert(tk.END, 
                        f"基因ID: {gene_id}\n"
                        f"匹配序列: {sequence}\n"
                        f"匹配位置: {match_info}\n\n"
                    )

    def reset_sequence(self):
        """重置当前序列"""
        if self.designer:
            self.designer.reset_sequence()
            # 重置后显示只有7-mer的序列
            self.update_interactive_display()

    def check_data_ready(self):
        """检查是否已准备好所需数据"""
        if not self.common_kmers or not self.off_targets:
            return False
        return True

    def reset_design_data(self):
        """重置设计相关数据"""
        self.common_kmers = {}
        self.off_targets = {}
        self.designer = None
        self.selected_kmer_entry.delete(0, tk.END)

    def use_existing_dict(self):
        """使用现有字典"""
        if not os.path.exists(OUTPUT_TEXT_FILE) or os.path.getsize(OUTPUT_TEXT_FILE) == 0:
            messagebox.showerror("错误", "未找到现有字典文件")
            return

        self.display_dictionary_content()
        messagebox.showinfo("成功", "已成功加载现有字典")

    def build_dict(self):
        """构建新字典"""
        if not self.file_path:
            messagebox.showerror("错误", "请先选择文件")
            return

        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, "分析中...\n")
        self.root.update()

        result = self.processor.build_dictionary(self.file_path)
        if result:
            self.result_text.insert(tk.END, 
                f"成功处理 {result['sequence_count']} 条序列\n"
                f"共找到 {result['kmer_count']} 个不同的7-mer\n"
            )
            self.display_dictionary_content()
        else:
            messagebox.showerror("错误", "构建字典失败")

    def filter_sequences(self):
        """筛选序列"""
        gene1_id = self.gene1_entry.get()
        gene2_id = self.gene2_entry.get()
        
        result = self.processor.process_gene_pair(gene1_id, gene2_id)
        if not result:
            messagebox.showerror("错误", "分析过程中出现错误")
            return

        # 保存结果到类属性中
        self.common_kmers, self.off_targets = result
        
        # 显示结果
        self.display_filter_results(self.common_kmers, self.off_targets)
        
        # 启用交互设计功能
        self.selected_kmer_entry.config(state='normal')

    def display_dictionary_content(self):
        """显示字典内容"""
        self.result_text.pack()
        self.result_text.delete(1.0, tk.END)
        with open(OUTPUT_TEXT_FILE, 'r') as f:
            self.result_text.insert(tk.END, f.read())

    def display_filter_results(self, common_kmers, off_target_genes):
        """显示筛选结果"""
        result_content = "共同的7-mers：\n"
        for kmer, details in common_kmers.items():
            result_content += f"\n{kmer}:\n"
            for detail in details:
                result_content += f"  {detail}\n"
            # 显示脱靶基因数量
            off_target_count = len(off_target_genes.get(kmer, {}).keys()) 
            result_content += f"  脱靶基因数量: {off_target_count}\n"

        result_content += "\n脱靶基因：\n"
        for kmer, gene_info in off_target_genes.items():
            result_content += f"\n{kmer}:\n"
            for gene_id, sequences in gene_info.items():  
                for seq_info in sequences:
                    sequence, position = seq_info
                    result_content += f"  {gene_id}: {sequence} ({position})\n"

        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, result_content)

    def show_file_selection(self):
        """显示文件选择相关组件"""
        self.entry.pack(pady=10)
        self.file_select_button.pack(pady=5)
        self.build_dict_button.pack(pady=5)
        self.result_text.pack(pady=10)

    def select_file(self):
        """选择文件"""
        self.file_path = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")]
        )
        self.entry_var.set(self.file_path)


