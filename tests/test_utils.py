import unittest
from core.utils import reverse_complement, highlight_sequence
from core.config import HIGHLIGHT_POSITIONS

class TestUtils(unittest.TestCase):
    def setUp(self):
        self.test_sequence = "ATGCTAGC"
        self.test_regions = [(2, 4), (6, 7)]

    def test_reverse_complement(self):
        """测试反向互补序列生成"""
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("AAATTT"), "AAATTT")
        self.assertEqual(reverse_complement(""), "")

    def test_highlight_sequence(self):
        """测试序列高亮标记"""
        result = highlight_sequence(self.test_sequence, self.test_regions)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]["sequence"], "TGC")
        self.assertEqual(result[1]["sequence"], "GC")

if __name__ == '__main__':
    unittest.main()