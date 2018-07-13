import os
import json
import unittest
import tempfile
import filecmp

from rna_blast_analyze.BR_core import convert_classes


class TestHTMLoutput(unittest.TestCase):
    def setUp(self):
        htmlfd, htmlo = tempfile.mkstemp()
        os.close(htmlfd)
        self.json_file = os.path.abspath(os.path.dirname(__file__) + '/test_data/RF00001_output.json')
        self.htmlo = htmlo

    def test_output(self):
        f = open(self.json_file, 'r')
        mydata = json.load(f)
        f.close()
        bb = convert_classes.blastsearchrecomputefromdict(mydata)
        bb.to_html(self.htmlo)
        try:
            self.assertTrue(
                filecmp.cmp(
                    self.htmlo,
                    os.path.abspath(
                        os.path.dirname(__file__) + '/test_data/RF00001_reference_output.html'
                    )
                )
            )
        finally:
            try:
                os.remove(self.htmlo)
            except:
                pass


if __name__ == '__main__':
    unittest.main()
