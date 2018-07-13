# functions for handling blast records, detecting and handling hits to same sequence segment


def merge_blast_hits(ext_seqs):
    # hits are grouped by sequence of origin
    # longer hits are first
    # or construct an dum sequence

    # change this to operate on procomputed sequence list, as manipulating sequences there is easier
    # sort new list by e-value

    # print('loading test data, remove before use')
    # for i, h in enumerate(ext_seqs):
    #     h.annotations['blast'][1].sbjct_start = 1500 + randint(0,50) + i
    #     h.annotations['blast'][1].sbjct_end = 2000 + randint(0,50) - i

    aligs = [i.annotations['blast'] for i in ext_seqs]

    alig_names = [i[0] for i in aligs]

    unique_orgns = set(alig_names)

    new_seqs_list = []

    for orgn_name in unique_orgns:
        if alig_names.count(orgn_name) == 1:

            q = [i for i, name in enumerate(alig_names) if name == orgn_name][0]
            new_seqs_list.append(ext_seqs[q])
            continue

        # the sequences should be iterated in order of decreasing length

        position_indexes = []
        overlaping_sequences = []

        same_source_seq_ind = [i for i, name in enumerate(alig_names) if name == orgn_name]

        hits = [aligs[ind] for ind in same_source_seq_ind]
        subset_ext_seqs = [ext_seqs[ind] for ind in same_source_seq_ind]

        sorted_hits = sorted(hits, key=_key_abs_seq_len_f, reverse=True)
        sorted_ext_seqs = sorted(subset_ext_seqs, key=_key_abs_seq_len_ext_f, reverse=True)

        for hit in sorted_hits:

            hit_ind = (hit[1].sbjct_start, hit[1].sbjct_end)
            overlap_with = check_ind_overlap(position_indexes, hit_ind)
            overlaping_sequences.append(overlap_with)
            position_indexes.append(hit_ind)

        # solve overlaping sequences
        non_overlap_rec = _solve_overlaping_sequeces(sorted_ext_seqs, overlaping_sequences)
        new_seqs_list += non_overlap_rec

    # sorted_by_eval = sorted(new_seqs_list, key=_key_eval_sort)
    sorted_by_original = reorder_seq_list_by_id(new_seqs_list, ext_seqs)

    if len(sorted_by_original) < len(ext_seqs):
        print('Total of {} Blast HSPs on input \n'
              '\t {} HSPs were found inside of longer HSP.\n'
              'They are removed from computation, as separate sequences, '
              'but will be included in output together with the longer HSP.\n'
              'Total sequences in computation: {}.'.format(
            len(ext_seqs),
            len(ext_seqs) - len(sorted_by_original),
            len(sorted_by_original)
        ))
        #todo solve overlaping hits output
        print('solve overlaping hits output')
    return sorted_by_original


def reorder_seq_list_by_id(target_list, sort_list):
    """
    function accepts list of SeqRecord objects and other list of either names or SeqRecord objects and reorders
     the first list by the second one
    :param target_list: list of SeqRecord objects to reorder
    :param sort_list: list of desired order
    :return:
    """
    # assert len(target_list) == len(sort_list)
    if hasattr(sort_list[0], 'id'):
        sort_by = [i.id for i in sort_list]
    else:
        sort_by = sort_list

    ns = []
    for sid in sort_by:
        for i, s in enumerate(target_list):
            if s.id == sid:
                ns.append(s)
                del target_list[i]
                break
    return ns


def _solve_overlaping_sequeces(ext_subset, overlaps):
    # overlaps is list of lists of overlaping sequences,
    # rows corespond to hits, nubers in the row coresponds to row number of sequence which the current one overlaps with

    new_list = []
    for ext, o in zip(ext_subset, overlaps):
        if not any(o):
            new_list.append(ext)
        else:
            longest_sequence_of_overlap = min(o)
            # it must already be in the new list
            new_list[longest_sequence_of_overlap].annotations['blast'].append(
                ext.annotations['blast'][1]
            )
            print('hits to same region found, merging hit to {}'.format(new_list[longest_sequence_of_overlap].id))

    return new_list


def _key_abs_seq_len_f(hit):
    # for sorting by blasti hit len on subject sequence
    # returns blast hit len on blast hsp
    return abs(hit[1].sbjct_start - hit[1].sbjct_end)


def _key_abs_seq_len_ext_f(ext_hit):
    # for sorting by blast hit len on subject sequence
    # returns blast hit len on extended sequence with blast hsp stored in annotations under 'blast' key
    return abs(ext_hit.annotations['blast'][1].sbjct_start - ext_hit.annotations['blast'][1].sbjct_end)


def _key_eval_sort(ext_hit):
    # for sorting by evalue
    return ext_hit.annotations['blast'][1].expect


def check_ind_overlap(ind, v):
    # only cares about if sequence is inside one of the indices in ind
    overlap = []
    for i, (s, e) in enumerate(ind):
        if s < e and v[0] < v[1]:
            if s <= v[0] and e >= v[1]:
                # inside overlap forward
                overlap.append(i)
            else:
                # hit not inside another
                pass

        elif s > e and v[0] > v[1]:
            if s >= v[0] and e <= v[1]:
                # inside overlap backward
                overlap.append(i)
            else:
                # hit not inside another
                pass

        else:
            # non-matching directions, checking for overlap not make sense
            pass
    return overlap
