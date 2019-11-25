import os
import json
import unittest
import tempfile
import hashlib

from rna_blast_analyze.BR_core import convert_classes
from rna_blast_analyze.BR_core.output.htmloutput import write_html_output


class TestHTMLoutput(unittest.TestCase):
    def setUp(self):
        htmlfd, htmlo = tempfile.mkstemp(prefix='rba_', suffix='_t11')
        os.close(htmlfd)
        self.json_file = os.path.abspath(os.path.dirname(__file__) + '/test_data/RF00001_output.json')
        self.htmlo = htmlo

    def test_output(self):
        f = open(self.json_file, 'r')
        mydata = json.load(f)
        f.close()
        bb = convert_classes.blastsearchrecomputefromdict(mydata)
        with open(self.htmlo, 'wb') as h:
            h.write(write_html_output(bb))
        try:
            with open(self.htmlo, 'rb') as f, open(
                os.path.abspath(
                    os.path.dirname(__file__) + '/test_data/RF00001_reference_output.html.md5'
                )
            ) as r:
                self.assertEqual(
                    hashlib.md5(f.read()).hexdigest(),
                    r.read()
                )
        finally:
            try:
                os.remove(self.htmlo)
            except:
                print('removing temporary test files failed')

    def test_output_with_sequence_fail(self):
        f = open(self.json_file, 'r')
        mydata = json.load(f)
        f.close()
        bb = convert_classes.blastsearchrecomputefromdict(mydata)
        hit = bb.hits.pop(1)
        hit.extension = None
        bb.hits_failed.append(hit)
        with open(self.htmlo, 'wb') as h:
            h.write(write_html_output(bb))

        try:
            with open(self.htmlo, 'rb') as f, open(
                os.path.abspath(
                    os.path.dirname(__file__) + '/test_data/RF00001_reference_missing_hit.html.md5'
                )
            ) as r:
                self.assertEqual(
                    hashlib.md5(f.read()).hexdigest(),
                    r.read()
                )
        finally:
            try:
                os.remove(self.htmlo)
            except:
                print('removing temporary test files failed')


if __name__ == '__main__':
    unittest.main()
