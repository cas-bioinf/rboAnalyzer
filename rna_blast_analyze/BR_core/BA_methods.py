import copy
import time
import pandas as pd
import logging
from copy import deepcopy

from rna_blast_analyze.BR_core.BA_support import Subsequences, get_hit_n, annotate_ambiguos_base

ml = logging.getLogger('rboAnalyzer')


class BlastSearchRecompute(object):
    """
    wrapper class for whole blast search
    """
    def __init__(self, args, query, iteration):
        self.hits = HitList()
        self.hits_failed = HitList()
        self.__all_hits = HitList()
        self.pandas = None
        self.query = query
        self.creation = time.time()
        self._runstat = 1
        self.args = args
        self.date_of_run = time.localtime()
        self.best_matching_model = None
        self.iteration = iteration
        self.multi_query = False
        self.msgs = []

    def copy_hits(self):
        hits = deepcopy(self.hits)
        for hit in hits:
            annotate_ambiguos_base(hit.extension)
        self.__all_hits = hits

    def get_all_hits(self):
        return deepcopy(self.__all_hits)

    def update_hit_stuctures(self):
        b_dict = {get_hit_n(h): h for h in self.hits}
        for h in self.hits:
            n = get_hit_n(h)
            if n in b_dict:
                h.extension.letter_annotations.update(b_dict[n].extension.letter_annotations)
                h.extension.annotations.update(b_dict[n].extension.annotations)

        for h in self.__all_hits:
            n = get_hit_n(h)
            if n in b_dict:
                h.extension.letter_annotations.update(b_dict[n].extension.letter_annotations)
                h.extension.annotations.update(b_dict[n].extension.annotations)

    def stop_timer(self):
        if self._runstat:
            self._runstat = 0
            self.total_time = time.time() - self.creation

    def to_pandas_dump(self, output_file):
        if self.pandas is None:
            self.export_pandas_results()
        self.pandas.to_pickle(output_file)

    def to_csv(self, output_file):
        if self.pandas is None:
            self.export_pandas_results()
        self.pandas.to_csv(output_file)

    def export_pandas_results(self):
        self.stop_timer()

        export_columns = [
            'blast_query',
            'subject',
            'bstart',
            'bend',
            'blast_bits',
            'best_sequence',
            'estart',
            'eend',
            'blast_eval',
            'query_start',
            'query_end',
            'b_e_start',
            'b_e_end',
            'locarna_score',
        ]

        # build export DataFrame
        # add prediction method names ad columns to final DataFrame (add lists to "data")
        #  each prediction method has its own columns
        #  unused columns are not listed
        data = dict()
        if isinstance(self.args.prediction_method, str):
            self.args.prediction_method = [self.args.prediction_method,]

        for col in export_columns + self.args.prediction_method:
            data[col] = []

        for i, hit in enumerate(self.hits):
            if hit.extension is None:
                ml.warning("Skipping exporting hit {}, which was not extended.".format(i))
                continue

            for k in data.keys():
                if k == 'blast_query':
                    # query name
                    data['blast_query'].append(self.query.id)
                    continue
                elif k == 'locarna_score':
                    data['locarna_score'].append(hit.extension.annotations['score'])
                    continue
                elif k == 'b_e_start':
                    data['b_e_start'].append(hit.source.annotations['extended_start'])
                    continue
                elif k == 'b_e_end':
                    data['b_e_end'].append(hit.source.annotations['extended_end'])
                    continue
                elif k == 'query_start':
                    data['query_start'].append(hit.source.annotations['blast'][1].query_start)
                    continue
                elif k == 'query_end':
                    data['query_end'].append(hit.source.annotations['blast'][1].query_end)
                    continue
                elif k == 'blast_eval':
                    data['blast_eval'].append(hit.source.annotations['blast'][1].expect)
                    continue
                elif k == 'blast_bits':
                    data['blast_bits'].append(hit.source.annotations['blast'][1].bits)
                    continue
                elif k == 'subject':
                    # subject sequence name
                    data['subject'].append(hit.source.id)
                    continue
                elif k == 'bstart':
                    # blast start
                    data['bstart'].append(hit.source.annotations['blast'][1].sbjct_start)
                    continue
                elif k == 'bend':
                    # blast end
                    data['bend'].append(hit.source.annotations['blast'][1].sbjct_end)
                    continue
                elif k == 'estart':
                    # extended start
                    data['estart'].append(hit.best_start)
                    continue
                elif k == 'eend':
                    # extended end
                    data['eend'].append(hit.best_end)
                    continue
                elif k == 'best_sequence':
                    # selected sequence
                    data['best_sequence'].append(str(hit.extension.seq))
                    continue
                elif hasattr(hit, k):
                    data[k].append(getattr(hit, k))
                    continue
                elif k in hit.extension.letter_annotations:
                    # add secondary structure
                    data[k].append(hit.extension.letter_annotations[k])
                else:
                    if k in self.args.prediction_method:
                        # key (prediction method name) is missing from letter annotations, but it was requested and
                        #  should be predicted. This means that something went wrong with the prediction and it was
                        #  logged with stdout. However, we need to include something to the table output
                        data[k].append('PREDICTION FAILED')
                    else:
                        raise Exception('key not found')
                pass
        self.pandas = pd.DataFrame(data)
        return self.pandas

    def write_results_fasta(self, out_file):
        self.stop_timer()
        with open(out_file, 'w') as f:
            for hit in self.hits:
                f.write('>{}\n{}\n'.format(
                    hit.extension.id,
                    str(hit.extension.seq))
                )

    def write_results_structures(self, out_file):
        with open(out_file, 'w') as f:
            for hit in self.hits:
                f.write('>{}\n{}\n'.format(
                    hit.extension.id,
                    str(hit.extension.seq))
                )
                for st_key in hit.extension.letter_annotations.keys():
                    f.write('{} {}\n'.format(hit.extension.letter_annotations[st_key], st_key))

    def res_2_record_list(self):
        return [hit.extension for hit in self.hits]

    def copy(self):
        new = BlastSearchRecompute(
            copy.deepcopy(self.args),
            self.query,
            self.iteration
        )
        new.hits = self.hits.copy()
        new.pandas = self.pandas
        new.multi_query = self.multi_query
        return new


class HitList(list):
    """
    Create list class where only subsequences objects can be stored
    overrides the append method
    """
    # def __init__(self, list):
    #     self.dpool_global = None
    #     self.dpool_raw = None
    def append(self, p_object):
        if not isinstance(p_object, Subsequences):
            raise Exception('passed object is not of expected class (Subsequences)')
        super(HitList, self).append(p_object)


def to_tab_delim_line_simple(input_args):
    A = vars(input_args)
    B = sorted(A.keys())
    line = '\t'.join([str(A[key]) for key in B])
    return line
