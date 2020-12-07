import unittest

from rna_blast_analyze.BR_core.predict_structures import encode_structure_unicode, decode_structure_unicode


class TestEncodeDecode(unittest.TestCase):
    def test_1(self):
        s = '.(((((((((((((.((((((((..((.(((...(((((((((((((...)))).....(((.((((....((((((((((..((((.(((.(((.....)))))).))))..)))))))))).)))).))).((((((....)))))))))))))))..))))).)))))))))))).)))))))))...'
        self.assertEqual(s, decode_structure_unicode(encode_structure_unicode(s)))


if __name__ == '__main__':
    unittest.main()