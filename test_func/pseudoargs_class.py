import json as json_package
import os
from argparse import Namespace


class Pseudoargs(Namespace):
    def __init__(
            self,
            query_in,
            blast_in,
            blast_db,
            b_type='guess',
            o_tbl=None,
            blast_regexp='(?<=\|)[A-Z0-9]*\.?\d*$',
            keep_all=False,
            max_shapes=5,
            subseq_step=3,
            subseq_window_locarna=30,
            subseq_window_muscle=30,
            subseq_window_mafft=50,
            subseq_window_simple_ext=10,
            reldev=(10, 20, 30),
            dump=None,
            pandas_dump=None,
            logfile=None,
            locarna_params='--struct-local 0 --sequ-local 0 --free-endgaps ++++',
            muscle_params='-distance1 kbit20_3 -distance2 pctidlog -objscore sp -weight1 threeway -weight2 none',
            mafft_params=' --quiet',
            mafft_mode='qinsi',
            mlocarna_params=None,
            locarna_anchor_length=7,
            repredict_file=None,
            pred_sim_threshold_percent=90,
            prediction_method=('alifold_unpaired_conserved_refold',),
            dill=None,
            pred_params=None,
            dev_pred=False,
            threads=1,
            ribosum=os.path.abspath(
                os.path.dirname(__file__) + '/../rna_blast_analyze/3rd_party_source/RSEARCH_matrices/RIBOSUM65.mat'
            ),
            structures2write=('best',),
            subseq_window=30,
            cm_file=None,
            use_rfam=False,
            config_file=None,
            html=None,
            json=None,
            csv=None,
            download_rfam=False,
            show_gene_browser=False,
            zip_json=False,
            **kwargs):
        super().__init__(**kwargs)
        self.blast_in = blast_in
        self.b_type = b_type
        self.o_tbl = o_tbl
        self.blast_db = blast_db
        self.blast_regexp = blast_regexp
        self.keep_all = keep_all
        self.max_shapes = max_shapes
        self.blast_query = query_in
        self.subseq_step = subseq_step
        self.subseq_window = subseq_window
        self.subseq_window_locarna = subseq_window_locarna
        self.subseq_window_mafft = subseq_window_mafft
        self.subseq_window_muscle = subseq_window_muscle
        self.subseq_window_simple_ext = subseq_window_simple_ext
        self.dump = dump
        self.pandas_dump = pandas_dump
        self.logfile = logfile
        self.locarna_params = locarna_params
        self.locarna_anchor_length = locarna_anchor_length
        self.repredict_file = repredict_file
        self.prediction_method = prediction_method
        self.dill = dill
        self.structures2write = structures2write
        self.ribosum = ribosum
        self.muscle_params = muscle_params
        self.mlocarna_params = mlocarna_params
        self.mafft_params = mafft_params
        self.mafft_mode = mafft_mode
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

        self.dev_pred = dev_pred

        if 0 < pred_sim_threshold_percent <= 100:
            self.pred_sim_threshold_percent = pred_sim_threshold_percent
        else:
            raise AttributeError('pred_sim_threshold_percent must be between 0 and 100')

        if isinstance(reldev, list):
            self.reldev = reldev
        else:
            self.reldev = list(reldev)

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