import os
import json
import logging
from jinja2 import Environment, FileSystemLoader, select_autoescape
from jinja2.exceptions import TemplateError
from time import strftime

from rna_blast_analyze.BR_core.BA_support import blasthsp2pre, run_rnaplot
from rna_blast_analyze.BR_core.config import CONFIG

from matplotlib import colors, cm

ml = logging.getLogger(__name__)


def rog_cmap():
    my_colors = [
        '#E24B2D',
        '#FFB916',
        '#7AD84B'
    ]
    cmap_rog = colors.LinearSegmentedColormap.from_list(
        name='rog',
        colors=[colors.hex2color(c) for c in my_colors],
        N=1024
    )
    return cmap_rog


def write_html_output(datain, template_path=''):
    ml.info("Writing HTML output.")
    # prepare data
    toprint = _prepare_body(datain)
    myfooter = _prepare_footer(datain)
    my_header = _prepare_header(datain)

    # init jinja2 rendering environment
    env = Environment(
        loader=FileSystemLoader(template_path),
        autoescape=select_autoescape(['html', 'xml'])
    )
    cwd = os.getcwd()
    try:
        os.chdir(CONFIG.html_template_dir)
        template = env.get_template('onehit.html')
        html_str = template.render(
            input_list=toprint,
            foo=myfooter,
            strftime=strftime,
            hea=my_header,
            show_gene_browser=datain.args.show_gene_browser,
            len=len,
        )
        return html_str
    except TemplateError:
        ml.error("Jinja rendering error. Please check if the template is available and correct.")
        raise
    except Exception:
        ml.error("Failed to render html.")
        raise
    finally:
        os.chdir(cwd)


def _prepare_header(data):
    return {
        'input': data.args.blast_in,
        'query': data.args.blast_query,
        'best_matching_model': data.best_matching_model,
    }


def _prep_hit_name(id, desc):
    if desc.startswith(id):
        return desc
    else:
        return id + " " + desc


def _prepare_body(data):
    rog = rog_cmap()
    norm = colors.Normalize(vmin=0, vmax=data.query.annotations['cmstat']['bit_sc']/2, clip=True)
    mm = cm.ScalarMappable(norm=norm, cmap=rog)

    jj = []
    for i, onehit in enumerate(data.hits):
        ext = onehit.extension
        h_bit_sc = ext.annotations['cmstat']['bit_sc']
        rr = dict()
        rr['source_seq_name'] = onehit.source.annotations['blast'][0]
        rr['seqname'] = ext.id
        rr['sequence'] = str(ext.seq)
        rr['formated_seq'] = ext.format('fasta')
        rr['rsearchbitscore'] = h_bit_sc
        rr['blast_hit_name'] = _prep_hit_name(onehit.source.annotations['blast'][0], onehit.source.description)
        rr['ext_start'] = onehit.best_start
        rr['ext_end'] = onehit.best_end
        rr['blast_text'] = blasthsp2pre(onehit.source.annotations['blast'][1])
        rr['pictures'] = _prepare_pictures(ext)
        rr['eval'] = onehit.source.annotations['blast'][1].expect
        rr['intid'] = str(i)
        rr['msgs'] = onehit.source.annotations['msgs']

        # estimate the homology
        q_sc = data.query.annotations['cmstat']['bit_sc']
        if h_bit_sc < 0:
            h_estimate = 'Not homologous'
        elif h_bit_sc/q_sc >= 0.5 and h_bit_sc >= 20:
            h_estimate = 'Homologous'
        else:
            h_estimate = 'Uncertain'

        rr['h_estimate'] = h_estimate
        rr['h_color'] = colors.rgb2hex(mm.to_rgba(h_bit_sc))

        # create seqviewurl here
        es = onehit.source.annotations['extended_start']
        ee = onehit.source.annotations['extended_end']
        if es > ee:
            es, ee = [ee, es]

        diff = 1000 + 2*len(ext)
        es -= diff
        ee += diff

        if es < 0:
            es = 1

        # if getattr(data.args, 'show_gene_browser', False):
        #     rr['draw_seqview'] = True
        # else:
        #     rr['draw_seqview'] = False

        rr['seqviewurl'] = ''.join([
            '?embeded=true',
            '&noviewheader=true',
            '&id={}'.format(onehit.source.annotations['blast'][0]),
            '&v={}:{}'.format(es, ee),
            '&appname=rboAnalyzer',
            '&mk={}:{}|BestMatch!'.format(onehit.best_start, onehit.best_end),
            '&multipanel=false',
            '&slim=true'
        ])
        rr['seqvid'] = 'seqv_{}'.format(i)

        jj.append(rr)
    return jj


def _prepare_pictures(sub):
    pictureslist = []
    for key in sub.letter_annotations.keys():
        np = dict()
        np['picname'] = key
        np['secondary_structure'] = sub.letter_annotations[key]

        picfile = run_rnaplot(
            seq=str(sub.seq),
            structure=sub.letter_annotations[key],
            format='svg'
        )
        with open(picfile) as f:
            np['pic'] = "data:image/svg+xml;utf8," + f.read()

        pictureslist.append(np)
        try:
            os.remove(picfile)
        except FileNotFoundError:
            print('cannot remove file: {}, file not found'.format(picfile))
        except OSError:
            print('cannot remove file: {}, file is directory'.format(picfile))

    return pictureslist


def _prepare_footer(data):
    # write command
    if hasattr(data.args, 'command') and data.args.command:
        command = ' '.join(data.args.command)
    else:
        command = 'run directly from python'

    # parameters
    params = [i for i in dir(data.args) if not i.startswith('__') and not callable(getattr(data.args, i))]
    p_text = []
    for arg in params:
        if arg == 'pred_params':
            js = json.dumps(
                getattr(data.args, arg),
                sort_keys=True,
                indent=4,
            )
            p_text.append((arg, js[1:-1]))
            continue
        p_text.append((arg, getattr(data.args, arg)))

    if not hasattr(data, 'date_of_run'):
        data.date_of_run = None

    prepared_footer_data = {
        'command': command,
        'parameters': p_text,
        'exec_date': data.date_of_run,
        'logdup': data.args.logmsgs
    }
    return prepared_footer_data
