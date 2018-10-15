import copy
import re
import time
import json
import pandas as pd
import gzip

from rna_blast_analyze.BR_core.BA_support import Subsequences
import rna_blast_analyze.BR_core.convert_classes
import rna_blast_analyze.BR_core.output.htmloutput


class BlastSearchRecompute(object):
    """
    wrapper class for whole blast search
    """
    def __init__(self):
        # self.dpool_global = None
        # self.dpool_raw = None
        self.hits = HitList()
        self.pandas = None
        self.query = None
        self.creation = time.time()
        self._runstat = 1
        self.args = None
        self.date_of_run = time.localtime()

    def stop_timer(self):
        if self._runstat:
            self._runstat = 0
            self.total_time = time.time() - self.creation

    def export_pandas_all(self):
        """
        wrap everything to the pandas datastructure, can latter use it for exports etc.
        :return:
        """
        self.stop_timer()
        for hit in self.hits:
            pass
        print('not implemented yet')

    def to_json(self, output_file, gzip_json=False):
        """
        write jsonnized object
        fully fledged reverse function is provided in
          rna_blast_analyze.BR_core.convert_classes.blastsearchrecomputefromdict
        :param output_file:
        :return:
        """

        j_obj = json.dumps(rna_blast_analyze.BR_core.convert_classes.blastsearchrecompute2dict(self))
        if gzip_json:
            with open(output_file + '.gz', 'wb') as ff:
                ff.write(
                    gzip.compress(j_obj.encode())
                )
        else:
            with open(output_file, 'w') as ff:
                ff.write(j_obj)


    def to_pandas_dump(self, output_file):
        if self.pandas is None:
            self.export_pandas_results()
        self.pandas.to_pickle(output_file)

    def to_csv(self, output_file):
        if self.pandas is None:
            self.export_pandas_results()
        self.pandas.to_csv(output_file)

    def export_pandas_results(self):
        # todo: add rsearch bits if computed
        #  add query name
        #  remove b_e_start b_e_end -> address this where neccessary
        #    if b_e_stat and b_e_end should be removed i need to compute the blast ext diff in eval_blast differently
        #    ex the values are used there and are not easily replaceable
        # todo revrite export so it exports data from one run to one file (ie multiple queries support)
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
            'bstrand',
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
        # todo: remove later - quick fix because some test code assigns prediction parameter for one method as str
        if isinstance(self.args.prediction_method, str):
            self.args.prediction_method = [self.args.prediction_method,]

        for col in export_columns + self.args.prediction_method:
            data[col] = []

        for hit in self.hits:
            for k in data.keys():
                if k == 'blast_query':
                    # query name
                    data['blast_query'].append(self.query.id)
                    continue
                elif k == 'locarna_score':
                    data['locarna_score'].append(hit.subs[hit.ret_keys[0]].annotations['score'])
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
                elif k == 'bstrand':
                    # blast strand
                    data['bstrand'].append(hit.source.annotations['blast'][1].strand)
                    continue
                elif k == 'best_sequence':
                    # selected sequence
                    data['best_sequence'].append(str(hit.subs[hit.ret_keys[0]].seq))
                    continue
                elif hasattr(hit, k):
                    data[k].append(getattr(hit, k))
                    continue
                elif k in hit.subs[hit.ret_keys[0]].letter_annotations:
                    # add secondary structure
                    data[k].append(hit.subs[hit.ret_keys[0]].letter_annotations[k])
                else:
                    if k in self.args.prediction_method:
                        # key (prediction method name) is missing from letter annotations, but it was requested and
                        #  should be predicted. This means that something went wrong with the prediction and it was
                        #  logged with stdout. However, we need to include something to tha table out
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
                f.write('>{}\n{}\n'.format(hit.subs[hit.ret_keys[0]].id,
                                           str(hit.subs[hit.ret_keys[0]].seq)))

    def write_results_structures(self, out_file):
        with open(out_file, 'w') as f:
            for hit in self.hits:
                f.write('>{}\n{}\n'.format(hit.subs[hit.ret_keys[0]].id,
                                           str(hit.subs[hit.ret_keys[0]].seq)))
                for st_key in hit.subs[hit.ret_keys[0]].letter_annotations.keys():
                    f.write('{} {}\n'.format(hit.subs[hit.ret_keys[0]].letter_annotations[st_key], st_key))

    def res_2_record_list(self):
        return [hit.subs[hit.ret_keys[0]] for hit in self.hits]

    def copy(self):
        new = BlastSearchRecompute()
        new.args = copy.deepcopy(self.args)
        new.hits = self.hits.copy()
        new.pandas = self.pandas
        new.query = self.query
        return new

    def extract_all_structures(self):
        all_structures = dict()
        for key in self.hits[0].annotaions['sss']:
            all_structures[key] = self.hits.extract_structures(key)
        return all_structures

    def to_html(self, outputfile):
        htmlstr = rna_blast_analyze.BR_core.output.htmloutput.write_html_output(self)
        with open(outputfile, 'w') as f:
            f.write(htmlstr)


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

    def extract_structures(self, str_key):
        """
        returns list of all structures from hits defined byt str_key
        :param str_key: key to structure for retrieving
        :return: list
        """
        outstr = []
        for seqr in self:
            outstr.append(seqr.letter_annotations[str_key])
        return outstr

    def iter_all_structures(self):
        """
        iterator for all structures defined in annotations['sss']
        structures will be returned as dict
        :return: dict
        """
        for seqr in self:
            yield {k: seqr.subs[seqr.ret_keys[0]].letter_annotations[k]
                   for k in seqr.subs[seqr.ret_keys[0]].annotations['sss']}


def to_tab_delim_line(input_args):
    A = vars(input_args)
    B = sorted(A.keys())
    if input_args.pred_params:
        pp = str(input_args.pred_params)
        l = [str(A[key]) for key in B]
        ol = []
        for p in input_args.pred_params:
            o = []
            flag = '|pred_params|' + re.sub(' ', '', str(p))
            for i in l:
                if i == pp:
                    i = re.sub(' ', '', str(p))
                if '.dump' in i:
                    spa = i.split('.')
                    i = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
                if '.pandas_dump' in i:
                    spa = i.split('.')
                    i = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
                if '.o_tbl' in i:
                    spa = i.split('.')
                    i = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
                if '.dill' in i:
                    spa = i.split('.')
                    i = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
                if '.pdf_out' in i:
                    spa = i.split('.')
                    i = '.'.join(spa[:-1]) + flag + '.' + spa[-1]
                o.append(i)
            ol.append('\t'.join(o))
        line = '\n'.join(ol)
    else:
        line = '\t'.join([str(i) for i in A.values()])
    return line


def to_tab_delim_line_simple(input_args):
    A = vars(input_args)
    B = sorted(A.keys())
    line = '\t'.join([str(A[key]) for key in B])
    return line


def to_tab_delim_header(input_args):
    A = vars(input_args)
    B = sorted(A.keys())
    return '\t'.join(B)


def add_loc_to_description(analyzed_hits):
    for hit in analyzed_hits.hits:
        d2a = '{}-{}'.format(hit.best_start, hit.best_end)
        # hit.source.description += d2a
        hit.subs[hit.ret_keys[0]].description += d2a
