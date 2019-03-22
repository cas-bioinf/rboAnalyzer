import json as json_package
import os
from argparse import Namespace
from rna_blast_analyze.BR_core.parse_accession import accession_regex


class Pseudoargs(Namespace):
    def __init__(
            self,
            blast_query,
            blast_in,
            blast_db,
            b_type='guess',
            blast_regexp=accession_regex,
            keep_all=False,
            subseq_window_locarna=30,
            subseq_window_simple_ext=10,
            dump=None,
            pandas_dump=None,
            logfile=None,
            locarna_params='--struct-local 0 --sequ-local 0 --free-endgaps ++++',
            locarna_anchor_length=7,
            repredict_file=None,
            prediction_method=('clustalo_alifold_unpaired_conserved_refold_rnafoldc',),
            pred_params=None,
            dev_pred=False,
            threads=1,
            cm_file=None,
            use_rfam=False,
            config_file=None,
            html=None,
            json=None,
            csv=None,
            download_rfam=False,
            show_gene_browser=False,
            zip_json=False,
            filter_by_eval=None,
            filter_by_bitscore=None,
            skipp_missing=False,
            enable_overwrite=False,
            db_type="blastdb",
            mode='locarna',
            verbose=0,
            **kwargs):
        super().__init__(**kwargs)
        self.blast_in = blast_in
        self.b_type = b_type
        self.blast_db = blast_db
        self.blast_regexp = blast_regexp
        self.keep_all = keep_all
        self.blast_query = blast_query
        self.subseq_window_locarna = subseq_window_locarna
        self.subseq_window_simple_ext = subseq_window_simple_ext
        self.dump = dump
        self.pandas_dump = pandas_dump
        self.logfile = logfile
        self.locarna_params = locarna_params
        self.locarna_anchor_length = locarna_anchor_length
        self.repredict_file = repredict_file
        self.prediction_method = prediction_method
        self.threads = threads
        self.cm_file = cm_file
        self.use_rfam = use_rfam
        self.config_file = config_file
        self.csv = csv
        self.json = json
        self.html = html
        self.download_rfam = download_rfam
        self.show_gene_browser = show_gene_browser
        self.zip_json = zip_json
        self.filter_by_eval = filter_by_eval
        self.filter_by_bitscore = filter_by_bitscore
        self.skip_missing = skipp_missing
        self.logmsgs = []
        self.dev_pred = dev_pred
        self.enable_overwrite = enable_overwrite
        self.db_type = db_type
        self.mode = mode
        self.verbose = verbose

        if self.filter_by_eval is not None and self.filter_by_bitscore is not None:
            raise AttributeError('filter_by_eval is not allowed with filter_by_bitscore')

        # manage default parameters
        default_file = os.path.abspath(
            os.path.dirname(__file__) + '/../rna_blast_analyze/BR_core/prediction_parameters.json'
        )
        if dev_pred:
            self.pred_params = pred_params
        else:
            params = dict()
            with open(default_file, 'r') as ff:
                params.update(json_package.load(ff))
            if pred_params:
                params.update(pred_params)
            self.pred_params = params