import unittest
from core.find_targets import find_common_kmers_and_off_targets

class TestFindTargets(unittest.TestCase):
    def setUp(self):
        # 创建测试用的基因-kmer字典
        self.test_gene_kmers = {
            'gene1': [('ATCGTAC', 'ATCGTACGTA', 1)],
            'gene2': [('ATCGTAC', 'TATCGTACGC', 2)],
            'gene3': [('ATCGTAC', 'GATCGTACA', 3)]
        }

    def test_find_common_kmers_and_off_targets(self):
        """测试共同kmers和脱靶基因查找"""
        common_kmers, off_targets = find_common_kmers_and_off_targets(
            self.test_gene_kmers, 'gene1', 'gene2'
        )
        
        # 验证共同kmers
        self.assertIn('ATCGTAC', common_kmers)
        
        # 验证脱靶基因
        self.assertIn('ATCGTAC', off_targets)
        self.assertIn('gene3', off_targets['ATCGTAC'])

if __name__ == '__main__':
    unittest.main()