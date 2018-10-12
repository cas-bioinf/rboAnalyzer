import os
import json
from rna_blast_analyze.BR_core import convert_classes
import tempfile
import hashlib


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

    target = os.path.abspath(os.path.dirname(__file__) + '/../RF00001_reference_output.html.md5')
    with open(html_file, 'rb') as f, open(target, 'w') as t:
        t.write(hashlib.md5(f.read()).hexdigest())

    os.remove(html_file)

if __name__ == '__main__':
    prepare_new_htmlout()
