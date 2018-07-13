import os
import json
from rna_blast_analyze.BR_core import convert_classes


def prepare_new_htmlout():
    json_file = os.path.abspath(os.path.dirname(__file__) + '/../RF00001_output.json')
    f = open(json_file, 'r')
    mydata = json.load(f)
    bb = convert_classes.blastsearchrecomputefromdict(mydata)
    bb.to_html(
        os.path.abspath(os.path.dirname(__file__) + '/../RF00001_reference_output.html')
    )


if __name__ == '__main__':
    prepare_new_htmlout()
