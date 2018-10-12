import os
import re
import subprocess
from tempfile import mkstemp
import logging

from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACUnambiguousRNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rna_blast_analyze.BR_core.parser_to_bio_blast import bread
from rna_blast_analyze.BR_core.config import CONFIG
from rna_blast_analyze.BR_core.retry_decorator import retry
from rna_blast_analyze.BR_core.fname import fname

ml = logging.getLogger(__name__)


@retry(ChildProcessError, tries=5)
def structures2shape_virtualbox(infile, out_file=None, shape_level=5, params=''):
    ml.info('Running structures2shape in virtualbox.')
    ml.debug(fname())
    ml.warning(
        'Warning, this is rather slow, especially for long structures it can take minutes to finish one structure,'
        ' if you need abstraction level 5, consider switching to db2shape python implementation'
        ' also for LSSU ends in segfault (from my experience)'
    )
    FNULL = open(os.devnull, 'w')

    if out_file:
        out = out_file
    else:
        fd, out = mkstemp()
        os.close(fd)

    virt_rapid_tempfile = virtualbox_tempfile()

    virtualbox_copy(infile, virt_rapid_tempfile)

    try:
        cmd = "sshpass -p {} ssh {} '{}RNAshapes -mode abstract --shapeLevel {} {} {}' > {}".format(
            CONFIG.SSH_PASS,
            CONFIG.SSH_USER,
            CONFIG.rnashapes_path,
            shape_level,
            params,
            virt_rapid_tempfile,
            out
        )
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = subprocess.call(cmd, shell=True)
        else:
            r = subprocess.call(
                cmd,
                stdout=FNULL,
                shell=True
            )
        if r:
            ml.error('structures2shapes failed')
            ml.error(cmd)
            raise ChildProcessError(r)
        out_shapes = _str2shape_parse(out)
    finally:
        FNULL.close()
        if not out_file:
            os.remove(out)

    virtualbox_delete(virt_rapid_tempfile)

    return out_shapes


@retry(ChildProcessError, tries=5)
def rapid_shapes_list_virtualbox(fasta_file, shape2use, rapidshapes_params='', out=None, shape_level=5):
    ml.info('Running rapidshapes in virtualbox.')
    ml.debug(fname())
    if shape2use == '':
        raise Exception('no shape found')
    if re.search('[^\[\] _]', shape2use):
        raise Exception('shapes can only be specified with "[", "]" and "_" characters')

    if out:
        outfile = out
    else:
        fd, outfile = mkstemp()
        os.close(fd)

    # make tempfile
    rapid_in_tempfile = virtualbox_tempfile()

    # copy fasta to mashine over scp
    virtualbox_copy(source=fasta_file, dest=rapid_in_tempfile)

    # run analysis
    FNULL = open(os.devnull, 'w')
    try:
        cmd = "sshpass -p {} ssh {} '{}RapidShapes {} --mode list --list ""{}"" --shapeLevel {} {}' > {}".format(
            CONFIG.SSH_PASS,
            CONFIG.SSH_USER,
            CONFIG.rnashapes_path,
            rapidshapes_params,
            shape2use,
            shape_level,
            rapid_in_tempfile,
            outfile
        )
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = subprocess.call(cmd, shell=True)
        else:
            r = subprocess.call(
                cmd,
                shell=True,
                stdout=FNULL,
                stderr=FNULL
            )

        if r:
            msgfail = 'Call to RapidShapes failed.'
            ml.error(msgfail)
            ml.error(cmd)
            ml.error('errorcode: {}'.format(r))
            raise ChildProcessError(msgfail)
    finally:
        FNULL.close()

    virtualbox_delete(rapid_in_tempfile)

    return outfile


def virtual_decorator(func):
    if func.__name__ == 'rapid_shapes_list' and CONFIG.USE_VIRTUAL:
        func = rapid_shapes_list_virtualbox
    elif func.__name__ == 'structures2shape' and CONFIG.USE_VIRTUAL:
        func = structures2shape_virtualbox
    else:
        pass

    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper


# predict suboptimal structures in some shape
@virtual_decorator
def rapid_shapes_list(fasta_file, shape2use, rapidshapes_params='', out=None, shape_level=5):
    """using rapid shapes list method
    really only for long rnas
    compilation takes forever - for 6S rna it is 6 sec per sequence

    # maybe the --binPath for the bellmansgap compiler is needed for successfull execution

    # consider capturing all output to suppres text output during prediction

    :param fasta_file: file with rnasequences
    :param shape2use: string of space separated shapes to use
    :param rapidshapes_params:
    :param out: output file path
    :param out: output file path
    :param shape_level:
    """

    # if isinstance(inshape, list) and all([isinstance(i, CastShape) for i in inshape]):
    #     shape2use = ' '.join([i.shape for i in inshape])
    # else:
    #     shape2use = inshape
    ml.info('Runing rapidshapes.')
    ml.debug(fname())
    if shape2use == '':
        raise Exception('no shape found')
    if re.search('[^\[\] _]', shape2use):
        raise Exception('shapes can only be specified with "[", "]" and "_" characters')

    if out:
        outfile = out
    else:
        fd, outfile = mkstemp()
        os.close(fd)

    FNULL = open(os.devnull, 'w')
    try:
        cmd = '{}RapidShapes {} --mode list --list "{}" --shapeLevel {} {} > {}'.format(
            CONFIG.rapidshapes_path,
            rapidshapes_params,
            shape2use,
            shape_level,
            fasta_file,
            outfile
        )
        ml.debug(cmd)
        if ml.getEffectiveLevel() == 10:
            r = subprocess.call(cmd, shell=True)
        else:
            r = subprocess.call(
                cmd,
                shell=True,
                stdout=FNULL,
                stderr=FNULL
            )

        if r:
            msgfail = 'Call to RapidShapes failed.'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)
    finally:
        FNULL.close()

    return outfile


@retry(ChildProcessError, tries=5)
def virtualbox_tempfile():
    ml.debug(fname())
    cmd = "sshpass -p {} ssh {} 'tfile=$(mktemp /tmp/virt_tmp.XXXXXXXX); echo $tfile'".format(
        CONFIG.SSH_PASS,
        CONFIG.SSH_USER,
        )
    try:
        ml.debug(cmd)
        report_str = subprocess.check_output(cmd, shell=True)
        return report_str.decode().strip()
    except:
        msgfail = 'make tempfile failed'
        ml.error(msgfail)
        ml.error(cmd)
        raise ChildProcessError(msgfail)


@retry(ChildProcessError, tries=5)
def virtualbox_copy(source, dest):
    ml.debug(fname())
    cmd = "sshpass -p {} scp {} {}:{}".format(
        CONFIG.SSH_PASS,
        source,
        CONFIG.SSH_USER,
        dest
    )
    try:
        ml.debug(cmd)
        r = subprocess.call(cmd, shell=True)
        if r:
            raise ChildProcessError('command: {} failed'.format(cmd))
    except:
        ml.error('Virtualbox copy failed.')
        ml.error(cmd)
        raise


@retry(ChildProcessError, tries=5)
def virtualbox_delete(file):
    ml.debug(fname())
    assert '-rf' not in file

    cmd = "sshpass -p {} ssh {} 'rm {}'".format(
        CONFIG.SSH_PASS,
        CONFIG.SSH_USER,
        file
    )
    ml.debug(cmd)
    r = subprocess.call(
        cmd,
        shell=True
    )
    if r:
        msgfail = 'Virtualbox delete failed.'
        ml.error(msgfail)
        ml.error(cmd)
        raise ChildProcessError(msgfail)
    return


def read_rapid_shapes_list_outfile(infile):
    """
    read rapidshapes --mode list file
    :param infile:
    :return:
    """
    seqs = []
    with open(infile, 'r') as txt:
        line = txt.readline()
        while line:
            if line[0] == '>':
                seq_name = line[1:].strip()

                line = txt.readline()
                sequence = line.split()[1]

                structures_info = []
                line = txt.readline()
                while not (line == '\n' or line == ''):
                    structures_info.append(line.split())
                    line = txt.readline()

                # build record
                new_rec = SeqRecord(Seq(sequence),
                                    id=seq_name)
                for si in structures_info:
                    new_rec.letter_annotations[si[3]] = si[1]
                    new_rec.annotations[si[3]] = (si[0], si[2])

                seqs.append(new_rec)

                del seq_name
                del sequence
                del structures_info
                del new_rec

            line = txt.readline()
    return seqs


class CastShape(object):
    """Shape object for RNAshapes -method cast parsing
    """
    def __init__(self):
        self.shape = ''
        self.seqs = []
        self.score = None
        self.mfe_ratio = None

    def shape_seq2str(self):
        """
        extract structures serving as templates
        :return:
        """
        if not self.seqs:
            raise NoSequencesError('No sequences found')
        qq = []
        for seq in self.seqs:
            if not isinstance(seq, SeqRecord):
                raise NoSequencesError('Assigned sequences in wrong format, expected SeqRecord from Biopython')
            qq.append({'str': seq.letter_annotations['ss0'], 'seq': str(seq.seq)})
        return qq


@virtual_decorator
def structures2shape(infile, out_file=None, shape_level=5, params=''):
    """
    convert structure to shape with RNAshapes
    file input must be fasta/like file but with structures instead of sequences
    structure cannot contain lonely pairs
    consider using db2shape function, which should give same result in most cases but is significantly faster
    """
    ml.info('Running Structures2shape')
    ml.debug(fname())
    ml.warning(
        'Warning, this is rather slow, especially for long structures it can take minutes to finish one structure,'
        ' if you need abstraction level 5, consider switching to db2shape python implementation'
        ' also for LSSU ends in segfault (from my experience)'
    )
    FNULL = open(os.devnull, 'w')

    if out_file:
        out = out_file
    else:
        fd, out = mkstemp()
        os.close(fd)

    try:
        cmd = '{}RNAshapes -mode abstract --shapeLevel {} {} {} > {}'.format(
            CONFIG.rnashapes_path,
            shape_level,
            params,
            infile,
            out
        )
        ml.debug(cmd)
        r = subprocess.call(
            cmd,
            stdout=FNULL,
            shell=True
        )
        if r:
            msgfail = 'Call to structures2shape failed'
            ml.error(msgfail)
            ml.error(cmd)
            raise ChildProcessError(msgfail)
        out_shapes = _str2shape_parse(out)
    finally:
        FNULL.close()
        if not out_file:
            os.remove(out)

    return out_shapes


def _str2shape_parse(fin):
    out = []
    with open(fin, 'r') as f:
        # each line is shape
        txt = f.readline()
        while txt:
            if re.search('[^\[\] _]', txt.strip()):
                raise Exception('shapes can only be specified with "[", "]" and "_" characters')
            r = CastShape()
            r.shape = txt.rstrip()
            out.append(r)
            txt = f.readline()
    return out


def cast_parse(file):
    """parse RNAshapes cast output
    output is grouped by distiguishible shapes and holds structure prediction for every input sequence
    returns list of cast_shapes objects
    """
    r_obj = []
    with open(file, 'r') as f:
        txt = bread(f)
        if 'No consensus shapes found' in txt:
            raise NoShapeFoundError('Common shape not found, try to increase relative or absolute energy deviation')
        while txt:
            if re.search('^\d+\)\s*(?=Shape)', txt):
                # enter a shape record
                curr_shape = CastShape()
                curr_shape.shape = _parse_shape(txt)
                curr_shape.score = _parse_score(txt)
                curr_shape.mfe_ratio = _parse_mfe_ratio(txt)

                curr_shape.seqs = _get_fasta_records(f)
                txt = bread(f)
                r_obj.append(curr_shape)
            else:
                raise UnboundLocalError(
                    'Unable to parse string "%s", expected str similar to "4)  '
                    'Shape: [][[][][][]][]  Score: -258.30  Ratio of MFE: 0.98"', txt
                )
    return r_obj


class Error(Exception):
    """Base class for exceptions here"""
    pass


class NoShapeFoundError(Error):
    def __init__(self, message):
        self.message = message


class TriedAllRelDevError(Error):
    def __init__(self, message):
        self.message = message


class NoSequencesError(Error):
    def __init__(self, message):
        self.message = message


def _get_fasta_records(f):
    rl = []
    txt = bread(f)
    while txt and txt[0] == '>':
        # fasta header
        [sid, sdescription] = _get_fasta_name(txt)
        txt = bread(f)
        # sequence record
        sr = re.search('[AUGCT]+', txt)
        # structure record
        txt = bread(f)
        st = txt[sr.start():sr.end()]

        # create sequence record
        record = SeqRecord(Seq(sr.group(), IUPACUnambiguousRNA),
                           id=sid, description=sdescription,
                           letter_annotations={'ss0': st}, annotations={'sss': ['ss0']})

        # todo Consider additional parsing for remaining fields of RNAshapes cast output.
        rl.append(record)
        txt = bread(f)
        if re.search('\d+\)', txt):
            # return one line
            f.seek(f.tell() - len(txt) - 2)
            break
    return rl


def _get_fasta_name(txt):
    """get id and definition from fasta header
    :return fasta_id (str) and fasta_description (str)"""
    # todo consider to rework this catching some more standard exceptions
    try:
        c = re.search('(?<=>) *[\S|.]+', txt)
        hit_id = c.group().lstrip()
        hit_def = txt[c.end() + 1:].lstrip()
    except:
        raise UnboundLocalError('Unable to parse fasta header from %s', txt)
    return hit_id, hit_def


def _parse_mfe_ratio(txt):
    """parse MFE ratio value from rnashapes cast entry header
    :return float"""
    # todo consider to rework this catching some more standard exceptions
    if ('Ratio' in txt or 'ratio' in txt) and 'of' in txt:
        try:
            rat = float(re.search('(?<=MFE:)\s*[01]\.\d+$', txt).group().lstrip().rstrip())
        except:
            raise UnboundLocalError('Unable to parse Ratio of "MFE" in string %s', txt)
    else:
        raise UnboundLocalError('Unable to parse "Ratio of" MFE in string %s', txt)
    return rat


def _parse_score(txt):
    """parse score value from rnashapes cast entry header
    :return float"""
    # todo consider to rework this catching some more standard exceptions
    try:
        score_val = float(re.search('(?<=Score:)\s*-?\d+\.?\d*\s*(?=Ratio)', txt).group().lstrip().rstrip())
    except:
        raise UnboundLocalError('Unable to parse a "Score" in string in %s', txt)
    return score_val


def run_rnashapes_cast(hits, reldev, keep_all=0):
    """run RNAshapes program mode cast
    input is """
    if len(hits) < 2:
        raise AssertionError('input to rnashapes must be at least 2 sequences')
    [fd, temp_name] = mkstemp()
    os.close(fd)
    err_name = None
    out_rnashapes_file = None

    # enable pop-ing th lovest value
    reldev.sort()
    reldev.reverse()
    try:
        seq_it = (seq for seq in hits)
        SeqIO.write(seq_it, temp_name, 'fasta')
        if len(reldev) == 0:
            raise TriedAllRelDevError('Tried all relative deviation levels from the list, '
                                      'please consider 3 parameter that might resolve this '
                                      '1) Provide more consise input fasta sequences. '
                                      '2) Decrease number of input fasta sequeces. '
                                      '3) Increase permisible relative deviation parameters (-reldev)')
        fd, out_rnashapes_file = mkstemp()
        err, err_name = mkstemp()
        call_shapes_cast(temp_name, reldev, fd, err)
        # call(['RNAshapes', '-mode','cast', '--relativeDeviation', str(reldev.pop()), temp_name],
        #      stdout=fd, stderr=err)

        # Result can be that no shapes are found
        # - so less sequences are prefered and ability for recomputation is needed
        with open(err_name, 'r') as e:
            err = e.read()
        if err:
            print(err)
            raise Exception(err)

        possible_shapes = cast_parse(out_rnashapes_file)
    except NoShapeFoundError:
        print('enter recrusion\n')
        # if not args.keep_all:
        #     os.remove(temp_name)
        #     os.remove(out_rnashapes_file)
        # run recursive
        possible_shapes = run_rnashapes_cast(hits, reldev)
    except:
        print('Something went wrong during execution of RNAshapes -mode cast, supported version is 3.3.0 (01.10.2015)')
        raise
    finally:
        if not keep_all:
            os.remove(temp_name)
            if err_name is not None:
                os.remove(err_name)
            if out_rnashapes_file is not None:
                os.remove(out_rnashapes_file)
    return possible_shapes


def call_shapes_cast(cast_infile, reldev, fd, err):
    """
    computes RNAshapes mode cast
    :param cast_infile: sequences infile
    :param reldev: permisible relative deviation
    :param fd: os level file descriptor for output file
    :param err: os level file descriptor for error file
    :return:
    """
    ml.info('Running RNAshapes cast')
    ml.debug(fname())
    cmd = [
        '{}RNAshapes'.format(CONFIG.rnashapes_path),
        '-mode',
        'cast',
        '--relativeDeviation',
        str(reldev.pop()),
        cast_infile
    ]
    r = subprocess.call(
        cmd,
        stdout=fd,
        stderr=err
    )
    if r:
        msgfail = 'call to RNAshapes mode cast failed'
        ml.error(msgfail)
        ml.error(cmd)
        raise ChildProcessError(msgfail)
    return


def _parse_shape(txt):
    """parse shape from rnashapes cast parse entry header
    :return string"""
    try:
        shape_string = re.search('(?<=Shape:)\s*[\[\]]+\s*(?=Score)', txt).group().lstrip().rstrip()
    except:
        raise UnboundLocalError('Unable to parse a "Shape" string in %s', txt)
    return shape_string
