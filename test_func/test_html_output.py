import os
import json
import unittest
import tempfile
import filecmp
import gzip

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
            # self.assertTrue(
            #     filecmp.cmp(
            #         self.htmlo + '.gz',
            #         os.path.abspath(
            #             os.path.dirname(__file__) + '/test_data/RF00001_reference_output.html.gz'
            #         )
            #     )
            # )
            with open(self.htmlo, 'rb') as f, gzip.open(
                os.path.abspath(
                    os.path.dirname(__file__) + '/test_data/RF00001_reference_output.html.gz'
                )
            ) as r:
                self.assertEqual(
                    f.read(),
                    r.read()
                )
        finally:
            try:
                os.remove(self.htmlo)
            except:
                print('removing temporary test files failed')


if __name__ == '__main__':
    unittest.main()
