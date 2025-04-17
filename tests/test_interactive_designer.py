import unittest
from core.interactive_designer import InteractiveDesigner
from core.config import MAX_MIRNA_LEN

class TestInteractiveDesigner(unittest.TestCase):
    def setUp(self):
        # 测试数据
        self.test_kmer = "ATCGTAC"
        self.test_common_kmers = {
            "ATCGTAC": ["gene1:ATCGTACGTA (1)", "gene2:TATCGTACGC (2)"]
        }
        self.test_off_targets = {
            "ATCGTAC": {
                "gene3": [("GATCGTACA", 3)]
            }
        }
        
        self.designer = InteractiveDesigner(
            self.test_kmer,
            self.test_common_kmers,
            self.test_off_targets
        )

    def test_initialization(self):
        """测试初始化"""
        self.assertEqual(self.designer.selected_kmer, self.test_kmer)
        self.assertEqual(self.designer.current_sequence, "")

    def test_add_base(self):
        """测试添加碱基"""
        self.designer.add_base("A")
        self.assertEqual(self.designer.current_sequence, "A")
        
        # 测试序列长度限制
        for _ in range(MAX_MIRNA_LEN):
            self.designer.add_base("A")
        self.assertLessEqual(len(self.designer.get_full_sequence()), MAX_MIRNA_LEN)

    def test_gc_content(self):
        """测试GC含量计算"""
        self.designer.add_base("G")
        self.designer.add_base("C")
        gc_content = self.designer.get_gc_content()
        # 计算期望的GC含量
        expected_gc = (4 / 9) * 100  # (2 + 2) / (7 + 2) * 100
        self.assertAlmostEqual(gc_content, expected_gc, places=2)

    def test_off_targets(self):
        """测试脱靶检测"""
        self.designer.add_base("A")
        count, genes = self.designer.get_off_targets()
        self.assertGreaterEqual(count, 0)
        self.assertIsInstance(genes, list)

    def test_reset_sequence(self):
        """测试序列重置"""
        self.designer.add_base("A")
        self.designer.reset_sequence()
        self.assertEqual(self.designer.current_sequence, "")

if __name__ == '__main__':
    unittest.main()