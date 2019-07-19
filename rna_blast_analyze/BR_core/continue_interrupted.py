import json
import logging
from rna_blast_analyze.BR_core.convert_classes import blastsearchrecomputefromdict, seqrecordfromdict
from rna_blast_analyze.BR_core.repredict_structures import wrapped_ending_with_prediction

ml = logging.getLogger('rboAnalyzer')


def continue_computation(intermediate_file):
    with open(intermediate_file, 'r') as f:
        data = json.load(f)
        analyzed_hits = blastsearchrecomputefromdict(data[0])
        new_structures = {k: [seqrecordfromdict(hit) for hit in v] for k, v in data[1].items()}
        exec_time = data[2]

    print('STATUS: Finishing computation for query: {}'.format(analyzed_hits.query.id))

    # get the cm model
    if analyzed_hits.args.cm_file:
        cm_file_rfam_user = analyzed_hits.args.cm_file
    else:
        cm_file_rfam_user = None

    wrapped_ending_with_prediction(
        args_inner=analyzed_hits.args,
        analyzed_hits=analyzed_hits,
        used_cm_file=cm_file_rfam_user,
        multi_query=analyzed_hits.multi_query,
        iteration=analyzed_hits.iteration,
        new_structures=new_structures,
        exec_time=exec_time
    )
    return analyzed_hits.iteration