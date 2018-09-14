import unittest
import os
from rna_blast_analyze.BR_core.parser_to_bio_blast import blast_parse_txt
from Bio.Blast import NCBIXML


class TestBlastParser(unittest.TestCase):
    def setUp(self):
        self.alig_attrs = [
            'hit_id',
            'length',
            'accession',
        ]

        self.hsp_attrs = [
            'query',
            'sbjct',
            'query_start',
            'query_end',
            'sbjct_start',
            'sbjct_end',
            # 'strand',
            # 'expect',
            # 'bits',
        ]

    def test_web_multi(self):
        blast_in = os.path.abspath(os.path.dirname(__file__) + '/test_data/web_multi_hit.txt')
        bixml = os.path.abspath(os.path.dirname(__file__) + '/test_data/web_multi_hit.xml')

        self._make_t(blast_in, bixml)

    def _make_t(self, file_txt, file_xml):
        with open(file_txt, 'r') as f, open(file_xml, 'r') as x:
            for btxt, bxml in zip(blast_parse_txt(f), NCBIXML.parse(x)):                                # queries
                for k in range(max((len(btxt.alignments), len(bxml.alignments)))):
                    self._check_attrs(btxt.alignments[k], bxml.alignments[k], self.alig_attrs)
                    for m in range(max((len(btxt.alignments[k].hsps), len(bxml.alignments[k].hsps)))):  # HSPs
                        self._check_attrs(
                            btxt.alignments[k].hsps[m],
                            bxml.alignments[k].hsps[m],
                            self.hsp_attrs
                        )

    def _check_attrs(self, i, j, attrlist):
        for attr in attrlist:
            try:
                if not hasattr(i, attr):
                    print("blast_parse dont have attribute: {}".format(attr))
                    raise AttributeError
                if not hasattr(j, attr):
                    print("blast_xml dont have attribute: {}".format(attr))
                    raise AttributeError
                assert getattr(i, attr) == getattr(j, attr)
            except AssertionError:
                if attr is not 'hit_id':                                         # hit_id is different
                    self.assertEqual(getattr(i, attr), getattr(j, attr))         # rerun assertion
                self._print_diff(i, j, attr)
            except AttributeError:
                if hasattr(i, attr):
                    print("blast_parse {}: {}".format(attr, getattr(i, attr)))
                if hasattr(j, attr):
                    print("blast_xml {}: {}".format(attr, getattr(j, attr)))
            except:
                raise

    @staticmethod
    def _print_diff(i, j, attr):
        if getattr(i, attr) != getattr(j, attr):
            print("blast_parse {}: {}".format(attr, getattr(i, attr)))
            print("blast_xml {}: {}".format(attr, getattr(j, attr)))


if __name__ == '__main__':
    unittest.main()
