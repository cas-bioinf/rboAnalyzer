# use bio blast classes but with my own parser only for nucleotide
import argparse
import copy
import logging
import math
import os
import re

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Blast import Record
from rna_blast_analyze.BR_core.fname import fname


ml = logging.getLogger(__name__)


def f_parser():
    """
    Input parser
    :return: args structure
    """
    parser = argparse.ArgumentParser(description='minimal parser for blast flat file')
    parser.add_argument('blast_in', help='blast file in')

    args = parser.parse_args()
    return args


# normaly, we have 2 methods, parse and read
# parse returns iterator on queries
# read assumes single query
# if multiple query on txt, we need to seek to EOF, several lines form the end parse the blast information, which is
# common to all query records

def blast_parse_txt(handle):
    """
    parser is not prepared to handle joined search outputs from different databases
        the parser do not check for this and will return all records as if they come from the database for last query
    it can handle search output for multiple queries to single database
    :param handle:
    :return:
    """
    ml.info(fname())
    lead_info = _read_leading_info(handle)
    trail_info = _read_trailing_info(handle)
    lead_info.update(trail_info)

    for one_query in _parse_blast_body(handle, lead_info):
        yield one_query


def _read_until_word_or_blank(txt, f, other_words=frozenset(('Database: ', 'Query=', 'Length=', 'RID: '))):
    line = []
    while txt:
        txt = f.readline()

        if not txt.strip():
            break
        if any(w in txt for w in other_words):
            break
        line.append(txt)
    return txt, ''.join(line)


def _read_leading_info(f):
    data = dict()

    txt = f.readline()
    if not txt:
        print('Your file {} appears to be empty'.format(f.name))
        raise IOError

    lt = f.tell()
    while not (txt[:6] == "Query="):
        lt = f.tell()
        if re.search('BLASTN \d+\.\d+\.\d+\+?', txt):
            m = re.search('BLASTN \d+\.\d+\.\d+\+?', txt)
            temp = m.group().split()
            data['application'] = temp[0]
            data['version'] = temp[1]

            txt = f.readline()

        elif re.search('(?<=RID: ).+', txt):
            m = re.search('(?<=RID: ).+', txt)
            data['rid'] = m.group()

            txt = f.readline()

        elif re.search('Reference:', txt):
            # reference spans several rows, usually blank line between
            txt, data['reference'] = _read_until_word_or_blank(txt, f)

        elif re.search('Database:', txt):
            # usualy two lines, end with "letters"
            txt, dbline = _read_until_word_or_blank(txt, f)
            al = dbline.split(";")
            tl = re.search('\d+(?:[\d,]*\d)', al[-1])
            if tl:
                data['database_letters'] = int(tl.group().replace(',', ''))
            ns = re.search('\d+ (?=sequences$)', al[-2])
            if ns:
                data['database_sequences'] = int(ns.group())
        else:
            txt = f.readline()
    f.seek(lt)
    return data


def _read_trailing_info(f):
    tail_txt = _get_trailing_section(f)
    return _parse_tail_section(tail_txt)


def _get_trailing_section(f):
    opos = f.tell()
    f.seek(0, os.SEEK_END)
    fl = f.tell()
    b = 1024

    def _get_ts(n):
        # f.seek(-b*n, os.SEEK_END)
        k = fl - b*n
        if k < 0:
            k = 0
        f.seek(k)
        return f.read(b*n)

    def _check_ts(part):
        # needed info
        #   Database:
        #   Lambda
        #   last subject (or) no hits found
        l = IUPACAmbiguousDNA.letters + 'U-'

        reg1 = re.compile("Database: ")
        reg2 = re.compile("Lambda\s+K\s+H\n?\s+\d+")
        reg3 = re.compile("(Sbjct\s+\d+\s+[" + l + "]+\s+\d+)|(\*{5}\s+No hits found\s+\*{5})")
        return all((
            re.search(reg1, part),
            re.search(reg2, part),
            re.search(reg3, part)
        ))

    block_i = 1
    section = _get_ts(block_i)

    while not _check_ts(section):
        block_i += 1
        section = _get_ts(block_i)

    f.seek(opos)

    return section


def _parse_tail_section(tail_text):
    # now section contains tail of the blast output
    # section contains Lambda, Database, cost matrix
    tt = tail_text.split("\n")

    data = dict()

    def _find_matches(pattern):
        ff = [m for m in [re.search(pattern, txt) for txt in tt] if m]
        assert len(ff) == 1
        return ff[0]

    def _find_index_match(pattern):
        return [i for i, m in enumerate([re.search(pattern, txt) for txt in tt]) if m]

    # matrix
    try:
        pat = re.compile("(?<=Matrix:)\s*\w+\s+matrix:?\s*-?\d+\s+-?\d")
        found = _find_matches(pat)
        ft = re.split(" |:", found.group().strip())

        assert len(ft) == 4
        data['matrix'] = ft[0]
        data['sc_match'] = int(ft[2])
        data['sc_mismatch'] = int(ft[3])
        del ft
        del pat
    except AssertionError:
        ml.debug('parsing data: matrix - failed')
    except:
        raise

    # gap penalties
    try:
        pat = re.compile("(?<=Gap Penalties:).*")
        found = _find_matches(pat)
        g1 = re.search("(?<=Existence:)\s*-?\d+", found.group())
        g2 = re.search("(?<=Extension:)\s*-?\d+", found.group())

        ft = [float(g.group().strip()) if g else float('NaN') for g in (g1, g2)]
        assert all([not math.isnan(i) for i in ft])
        data['gap_penalties'] = ft
        del ft
        del pat
        del g1
        del g2
    except AssertionError:
        ml.debug('parsing data: gap_penalties - failed')
    except:
        raise

    # lambda k, h
    try:
        p1 = re.compile("Lambda\s+K\s+H")
        p2 = re.compile("Gapped")
        p1m = _find_index_match(p1)
        p2m = _find_index_match(p2)

        # non-gapped and gapped version
        # gapped is preceded by "Gapped" in prev line
        assert len(p1m) == 2
        p1r = [float(i) for i in tt[p1m[0]+1].strip().split()]
        assert len(p1r) == 3
        assert p2m[0] + 1 == p1m[1]

        p2r = [float(i) for i in tt[p1m[1]+1].strip().split()]
        assert len(p2r) == 3
        data['ka_params'] = tuple(p1r)
        data['ka_params_gap'] = tuple(p2r)
        del p1, p2
        del p1m, p2m
        del p1r, p2r

    except AssertionError:
        ml.debug('parsing data: ka_params - failed')
    except:
        raise

    # battery format xx: num
    bat = {
        "num_sequences": "Number of Sequences:",  # duplicate
        "num_hits": "Number of Hits to DB:",
        "num_good_extends": "Number of extensions:",
        "num_letters_in_database": "Length of database:",
        "num_sequences_in_database": "Number of sequences in database:",  # duplicate
        "effective_search_space_used": "Effective search space used:",
        "hsps_gapped": "Number of HSP's gapped:",
        "effective_search_space": "Effective search space:",
        "effective_query_length": "Effective length of query:",
        "effective_database_length": "Effective length of database:"
    }

    for key in bat.keys():
        try:
            pat = re.compile("(?<=" + bat[key] + ")\s*\d+")
            m = _find_matches(pat)
            ft = m.group().strip()
            data[key] = int(ft)
            del pat
            del m
            del ft
        except AssertionError:
            ml.debug("parsing data: {} - failed".format(key))
        except:
            raise
    del bat

    # better then
    bat = {
        "num_seqs_better_e": "Number of sequences better than",
        "hsps_no_gap": "Number of HSP's better than"
    }
    for key in bat.keys():
        try:
            # capture all, then do the easy part, then match the posible eval
            pat = re.compile("(?<=" + bat[key] + ").*")
            m = _find_matches(pat)

            # capture count
            ft = m.group().split()
            assert len(ft) > 1
            data[key] = int(ft[-1].split(':')[-1].strip())  # try split with ":" it may not be necessary

        except AssertionError:
            ml.debug("parsing data: {} - failed".format(key))
        except:
            raise

    return data


def update_blast_record(br, commonparseddata):
    """
    :param br: Record
    :param commonparseddata: dict
    :return:
    """
    for key in commonparseddata.keys():
        setattr(br, key, commonparseddata[key])


def _parse_blast_body(f, common_info):
    """
    run on open handle to file
    for each query return list of blast records

    web format differ from standalone export format
    standallone 2.2.8+
        output is segmented by query, each query is ended with KH report
        contain all queries, if no hit for given query: "***** No hits found *****"

    web:
        output is segmented by queries, new query is signaled only by "Query="
        if no hits are found, the query is not part of the report and not signaled in any way that that it was in the
         input
        if no hits are returned blast will not allow txt file download
    """

    # start the parse loop
    txt = bread(f)

    record_holder = Record.Blast()
    update_blast_record(record_holder, common_info)
    do_query_name = True
    do_query_length = True
    R = copy.deepcopy(record_holder)
    while txt != '':
        # parse query name
        if do_query_name and re.search('^Query=', txt):
            # query_name = re.search('(?<=Query=)\ *\S+', txt).group().lstrip()
            # it is possible, that no query name is provided (web)
            R.query = re.search('(?<=Query=)\s*\S*', txt).group().lstrip()
            do_query_name = False
            txt = bread(f)
            continue

        # parse query Lenght
        if do_query_length and re.search('^Length=', txt):
            R.query_length = int(re.search('(?<=Length=)\ *\d+', txt).group().lstrip())
            txt = bread(f)
            do_query_length = False
            continue

        # skip reference, database and table
        # skip table
        # parse hits to sequence
        if txt[:10] == 'ALIGNMENTS' or txt[0] == '>':
            # enter the alignments loop
            if txt[0] == '>':
                # seek before the line as the subsequent function expects it
                f.seek(f.tell() - len(txt) - 2)
            [R, txt] = read_aligns(f, R)
            # R.query_length = query_length
            # R.query = query_name
            # append blast record object to the temporary wrapper
            yield R
            R = copy.deepcopy(record_holder)
            do_query_name = True
            do_query_length = True
        elif re.search('\*{5}\sNo hits found\s\*{5}', txt) or re.search("^Query=", txt):
            # new record
            yield R
            R = copy.deepcopy(record_holder)
            do_query_name = True
            do_query_length = True
            txt = bread(f)
        else:
            # get the next line
            txt = bread(f)


def read_aligns(f, R):
    """
    read alignments
    parse all aligns to end of blast report or to next query sequence

    issue:
     now, genebank switched to using accession instead of gis
     and in blast there are accession.version identifier (and sometimes something different)

    """
    # parse_regexp = '[ACTGUactgu\-]+'
    parse_regexp_q = '[ACGTUKSYMWRBDHVNacgtuksymwrbdhvn\-]+'
    eor = True
    txt = bread(f)  # get line after alignments
    # if txt[0] == '>':  # enter a record
    #     # create an alignment object if we are sure that there is one
    #     # start of parsing hits to single source sequence
    while txt and eor:
        if txt[:6] == 'Query=':
            break
        c = re.search('(?<=>) *[\S|.]+', txt)
        def_rem = bread(f)
        while not re.search('Length=', def_rem):
            txt += def_rem
            def_rem = bread(f)
        # get new alig object
        curr_alig = Record.Alignment()

        curr_alig.hit_id = c.group()
        curr_alig.hit_def = txt[c.end() + 1:]
        curr_alig.length = int(re.search('(?<=Length=)\d+', def_rem).group())

        # draw next line
        txt = bread(f)

        while txt and eor:
            # loop the hsps of the record
            if txt[0] == '>' or txt[:6] == 'Query=':
                # R.alignments.append(curr_alig)
                break

            curr_hsp = Record.HSP()

            while txt and eor:
                # parse scores
                if txt[:8] == " Score =":
                    # changing "score" to bits
                    # because the score is in "bits" units
                    curr_hsp.bits = float(
                        re.search('(?<=Score =) *\d+\.?\d*', txt).group(0).lstrip()
                    )
                    # if re.search('bits \(', txt):

                    # changing "bits" to score, because the number in brackets is match score (in nucleotide blast
                    #  same as number of identities)
                    curr_hsp.score = int(
                        re.search('(?<=bits \() *\d+(?=\))', txt).group(0).lstrip()
                    )
                    # the regexp for scientific format from here:
                    # http://stackoverflow.com/questions/18152597/extract-scientific-number-from-string
                    curr_hsp.expect = float(
                        re.search('(?<=Expect =)-? *[0-9]+\.?[0-9]*(?:[Ee] *-? *[0-9]+)?', txt).group(0).lstrip()
                    )
                elif txt[:13] == " Identities =":
                    curr_hsp.identities = tuple(
                        [
                            int(i) for i in re.search(
                                '(?<=Identities =) *\d+/\d+(?= *\()', txt
                            ).group(0).lstrip().split('/')
                        ]
                    )
                    curr_hsp.gaps = tuple(
                        [
                            int(i) for i in re.search(
                                '(?<=Gaps =) *\d+/\d+(?= *\()', txt
                            ).group(0).lstrip().split('/')
                        ]
                    )
                    curr_hsp.align_length = curr_hsp.identities[1]
                elif txt[:8] == " Strand=":
                    curr_hsp.strand = tuple(
                        re.search(
                            '(?<=Strand=)\S+', txt).group(0).rstrip().split('/')
                    )
                elif txt[:10] == " Features ":
                    # must be enabled even if output isn't used, moves the pointer pas the features field
                    parsed_features = _parse_features(f, txt)
                    curr_hsp.features = parsed_features

                elif txt[:5] == "Query":
                    break

                else:
                    ml.debug("line ignored pos {} - {}".format(f.tell(), txt))

                # read next line
                txt = bread(f)

            # start of alignment (hsps) parse block
            # hsps starts with "Score ="
            qseq = ''
            mid = ''
            sseq = ''
            count = 1
            # while not re.search('Score =', txt):
            while txt[:8] != " Score =" or txt[:10] != " Features ":
                # Features in this part
                # Features flanking this part

                # parsing individual hps alignments
                # alignment may:
                #   1) continue with another line triple signaled by Query
                #   2) break
                #       a) continue with next hsps
                #       b) continue with next alignment (different source sequence (organism))
                #       c) be an EOF (this is not signaled) - tail was already parsed

                # read alignment
                # read alignment by triples
                if not txt[:5] == 'Query':
                    # this is eof (may not be)
                    while txt != '':
                        # if 'Database' in txt or 'database' in txt:
                        #     pass
                        if (txt[:6] == 'Query=') or (txt[0] == ">"):
                            # raise AssertionError('parser failed - unparsed record at %s', txt)
                            break
                        txt = bread(f)
                    # eor = False
                    break

                # get query start only at first instance
                if count == 1:
                    query_start = int(re.search('\d+', txt[5:]).group())

                # match query end at each instance
                query_end = int(re.search('\d+$', txt).group())

                # allow lowercase masking sequences
                # add 5 to prevent matching in Query
                q_info = re.search(parse_regexp_q, txt[5:])
                qseq += q_info.group()

                # get middle line
                txt = bread(f)
                mid += txt[q_info.start() + 5:q_info.end() + 5]

                txt = bread(f)
                sseq += txt[q_info.start() + 5:q_info.end() + 5]

                if count == 1:
                    subject_start = int(re.search('\d+', txt[5:]).group())

                subject_end = int(re.search('\d+$', txt).group())

                # go next
                txt = bread(f)
                count += 1
                if (txt[0] == '>') or (txt[:8] == " Score =") or (txt[:6] == 'Query=') or txt[:10] == " Features ":
                    break

            # if end of iteration save current hsp and go for next one
            curr_hsp.query = qseq
            curr_hsp.match = mid
            curr_hsp.sbjct = sseq
            curr_hsp.query_start = query_start
            curr_hsp.query_end = query_end
            curr_hsp.sbjct_start = subject_start
            curr_hsp.sbjct_end = subject_end

            assert len(qseq) == len(mid) == len(sseq)

            # append hsps to current alignment
            curr_alig.hsps.append(curr_hsp)

        R.alignments.append(curr_alig)

    if not eor and len(R.alignments) == 0:
        R.alignments.append(curr_alig)
    if len(R.alignments) == 0:
        print('no hits parsed')
    return R, txt


def _parse_features(f, txt):
    """
    there is two kinds of features
    1) Features in this part of subject sequence:
    2) Features flanking this part of subject sequence:

    in both cases, they may be multiline

    :param file_handle:
    :param txt: str
    :return:
    """

    pf = None
    if txt[:43] == " Features in this part of subject sequence:":
        pf = _features_inside(f, f.readline())
    elif txt[:49] == " Features flanking this part of subject sequence:":
        pf = _features_flanking(f, f.readline())
    else:
        # just find the start of expected line
        t = f.tell()
        while not (txt[:5] == "Query"):
            t = f.tell()
            txt = f.readline()
        # seek to before query was encountered
        f.seek(t)
    return pf


def _features_inside(f, txt):
    feature = []
    while txt[:3] == "   ":
        feature.append(txt)
        txt = f.readline()

    return "".join(feature)


def _parse_features_flank(txt):
    return int(txt.strip().split()[0]), "side".join(txt.split("side:")[1:])


def _features_flanking_strict(f, txt):
    flanks = [None, None]
    while txt[:3] == "   ":
        if re.search("at 5' side:", txt):
            flanks[0] = _parse_features_flank(txt)

        elif re.search("at 3' side:", txt):
            flanks[1] = _parse_features_flank(txt)
        else:
            ml.debug("flanking feature parsing failed for line: {}".format(txt))

        txt = f.readline()
    return flanks


def _features_flanking(f, txt):
    flanks = []
    while txt[:3] == "   ":
        if re.search("at 5' side:", txt):
            flanks.append(txt)
        elif re.search("at 3' side:", txt):
            flanks.append(txt)
        else:
            ml.debug("flanking feature parsing failed for line: {}".format(txt))

        txt = f.readline()
    return flanks


def bread(f):
    """read one line and skip blank lines
    return empty line when EOF"""

    t = f.readline()
    while t == '\n':
        t = f.readline()
    return t
