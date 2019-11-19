import os
import json
from rna_blast_analyze.BR_core import convert_classes
from rna_blast_analyze.BR_core.output.htmloutput import write_html_output
import tempfile
import hashlib


def prepare_new_htmlout():
    json_file = os.path.abspath(os.path.dirname(__file__) + '/../RF00001_output.json')
    f = open(json_file, 'r')
    mydata = json.load(f)
    bb = convert_classes.blastsearchrecomputefromdict(mydata)
    fd, html_file = tempfile.mkstemp(prefix='rba_', suffix='_t30')
    os.close(fd)

    with open(html_file, 'wb') as h:
        h.write(write_html_output(bb))

    target = os.path.abspath(os.path.dirname(__file__) + '/../RF00001_reference_output.html.md5')
    with open(html_file, 'rb') as f, open(target, 'w') as t:
        t.write(hashlib.md5(f.read()).hexdigest())

    os.remove(html_file)

    bb.hits[1].extension = None
    with open(html_file, 'w') as h:
        h.write(write_html_output(bb))

    target = os.path.abspath(os.path.dirname(__file__) + '/../RF00001_reference_missing_hit.html.md5')
    with open(html_file, 'rb') as f, open(target, 'w') as t:
        t.write(hashlib.md5(f.read()).hexdigest())
    os.remove(html_file)


if __name__ == '__main__':
    prepare_new_htmlout()
