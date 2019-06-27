import unittest
import operator
from rna_blast_analyze.BR_core.BA_verify import version_check


class TestBlastParser(unittest.TestCase):
    def test_less(self):
        self.assertTrue(
            version_check([2, 1, 0], [2, 0, 0], '', '')
        )

    def test_greater(self):
        self.assertTrue(
            version_check([2, 0, 0], [2, 1, 0], '', '', op=operator.le)
        )

    def test_equal(self):
        self.assertTrue(
            version_check([1, 1, 1], [1, 1, 1], '', '')
        )

    def test_equal2(self):
        self.assertTrue(
            version_check([1, 1, 1], [1, 1, 1], '', '', op=operator.le)
        )