import os
import json
from rna_blast_analyze.BR_core import convert_classes
import gzip
import tempfile


def prepare_new_htmlout():
    json_file = os.path.abspath(os.path.dirname(__file__) + '/../RF00001_output.json')
    f = open(json_file, 'r')
    mydata = json.load(f)
    bb = convert_classes.blastsearchrecomputefromdict(mydata)
    fd, html_file = tempfile.mkstemp()
    os.close(fd)
    bb.to_html(
        html_file
    )
    target = os.path.abspath(os.path.dirname(__file__) + '/../RF00001_reference_output.html.gz')
    with open(html_file, 'rb') as f, gzip.open(target, 'wb') as t:
        t.write(f.read())

    os.remove(html_file)

if __name__ == '__main__':
    prepare_new_htmlout()
