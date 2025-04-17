import unittest
from unittest.mock import patch, MagicMock
from core.processor import KmerProcessor
from core.config import OUTPUT_TEXT_FILE, OUTPUT_RESULT

class TestKmerProcessor(unittest.TestCase):
    def setUp(self):
        """设置测试环境"""
        self.processor = KmerProcessor()
        self.test_kmers = {
            'ATCGTAC': ['gene1:ATCGTACGTA (1)', 'gene2:TATCGTACGC (2)']
        }
        self.test_genes = {
            'gene1': [('ATCGTAC', 'ATCGTACGTA', 1)],
            'gene2': [('ATCGTAC', 'TATCGTACGC', 2)]
        }

    @patch('core.processor.build_kmer_dictionary')
    @patch('core.processor.save_to_text')
    def test_build_dictionary(self, mock_save, mock_build):
        """测试字典构建"""
        # 设置 mock
        mock_build.return_value = (self.test_kmers, self.test_genes)
        mock_save.return_value = True
        
        # 执行测试
        result = self.processor.build_dictionary("test.fasta")
        
        # 验证结果
        self.assertIsNotNone(result)
        self.assertEqual(result['sequence_count'], len(self.test_genes))
        self.assertEqual(result['kmer_count'], len(self.test_kmers))
        
        # 验证调用
        mock_build.assert_called_once_with("test.fasta")
        mock_save.assert_called_once_with(self.test_kmers, self.test_genes, OUTPUT_TEXT_FILE)

    @patch('core.processor.load_from_text')
    @patch('core.processor.find_common_kmers_and_off_targets')
    def test_process_gene_pair(self, mock_find, mock_load):
        """测试基因对处理"""
        # 设置 mock
        mock_load.return_value = (self.test_kmers, self.test_genes)
        mock_find.return_value = (
            {'ATCGTAC': ['gene1:seq1', 'gene2:seq2']},
            {'ATCGTAC': {'gene3': [('seq3', 3)]}}
        )
        
        # 执行测试
        result = self.processor.process_gene_pair('gene1', 'gene2')
        
        # 验证结果
        self.assertIsNotNone(result)
        common_kmers, off_targets = result
        self.assertTrue(isinstance(common_kmers, dict))
        self.assertTrue(isinstance(off_targets, dict))
        
        # 验证调用
        mock_load.assert_called_once_with(OUTPUT_TEXT_FILE)
        mock_find.assert_called_once_with(self.test_genes, 'gene1', 'gene2')

    @patch('core.processor.load_from_text')
    def test_process_gene_pair_with_invalid_data(self, mock_load):
        """测试无效数据情况"""
        mock_load.return_value = (None, None)
        result = self.processor.process_gene_pair('gene1', 'gene2')
        self.assertIsNone(result)

if __name__ == '__main__':
    unittest.main()