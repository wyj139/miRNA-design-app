import unittest
from unittest.mock import patch, MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from core.build_kmer_dict import build_kmer_dictionary

class TestBuildKmerDict(unittest.TestCase):
    def setUp(self):
        """设置测试数据"""
        # 创建测试序列记录
        self.test_sequences = [
            SeqRecord(
                Seq("ATGCTAGCTAGCTGATGC"),
                id="test_gene_1",
                description="test_gene_1"
            ),
            SeqRecord(
                Seq("CGTAGCTAGCTAGCTGAT"),
                id="test_gene_2",
                description="test_gene_2"
            )
        ]

    @patch('core.build_kmer_dict.read_fasta_file')
    def test_build_kmer_dictionary(self, mock_read_fasta):
        """测试字典构建功能"""
        # 设置 mock 返回值
        mock_read_fasta.return_value = self.test_sequences
        
        # 调用测试函数
        kmer_indices, gene_kmers = build_kmer_dictionary("test.fasta")
        
        # 验证结果
        self.assertIsNotNone(kmer_indices)
        self.assertIsNotNone(gene_kmers)
        self.assertTrue(isinstance(kmer_indices, dict))
        self.assertTrue(isinstance(gene_kmers, dict))
        
        # 验证内容
        self.assertGreater(len(kmer_indices), 0)
        self.assertGreater(len(gene_kmers), 0)
        
        # 验证数据格式
        for kmer, indices in kmer_indices.items():
            self.assertEqual(len(kmer), 7)  # 验证 kmer 长度
            self.assertTrue(all(isinstance(idx, str) for idx in indices))
            
        for gene_id, kmers in gene_kmers.items():
            self.assertTrue(all(len(k[0]) == 7 for k in kmers))  # 验证 kmer 长度
            self.assertTrue(all(isinstance(k[1], str) for k in kmers))  # 验证序列
            self.assertTrue(all(isinstance(k[2], int) for k in kmers))  # 验证位置

    def test_build_kmer_dictionary_with_invalid_file(self):
        """测试无效文件情况"""
        result = build_kmer_dictionary("invalid_file.fasta")
        self.assertEqual(result, (None, None))

    @patch('core.build_kmer_dict.read_fasta_file')
    def test_build_kmer_dictionary_with_progress(self, mock_read_fasta):
        """测试进度回调功能"""
        mock_read_fasta.return_value = self.test_sequences
        progress_values = []
        
        def progress_callback(value):
            progress_values.append(value)
            
        build_kmer_dictionary("test.fasta", progress_callback)
        
        self.assertTrue(len(progress_values) > 0)
        self.assertEqual(progress_values[-1], 100.0)

if __name__ == '__main__':
    unittest.main()