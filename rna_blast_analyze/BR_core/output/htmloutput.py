import os
from jinja2 import Environment, FileSystemLoader, select_autoescape
from time import strftime

from rna_blast_analyze.BR_core.BA_support import blasthsp2pre, run_rnaplot
from rna_blast_analyze.BR_core.config import CONFIG


def write_html_output(datain, template_path=''):
    # prepare data
    toprint = _prepare_body(datain)
    myfooter = _prepare_footer(datain)

    # init jinja2 rendering enviroment
    env = Environment(
        loader=FileSystemLoader(template_path),
        autoescape=select_autoescape(['html', 'xml'])
    )
    cwd = os.getcwd()
    try:
        os.chdir(CONFIG.html_template_dir)
        template = env.get_template('onehit.html')
        html_str = template.render(input_list=toprint, foo=myfooter, strftime=strftime)
    finally:
        os.chdir(cwd)
    return html_str


def _prepare_body(data):
    jj = []
    for i, onehit in enumerate(data.hits):
        ext = onehit.subs[onehit.ret_keys[0]]
        rr = dict()
        rr['source_seq_name'] = onehit.source.annotations['blast'][0]
        rr['seqname'] = ext.id
        rr['sequence'] = str(ext.seq)
        rr['formated_seq'] = ext.format('fasta')
        rr['rsearchbitscore'] = ext.annotations['cmstat']['bit_sc']
        rr['blast_hit_name'] = onehit.source.annotations['blast'][0] + ' ' + ' '.join(onehit.source.description.split()[1:])
        rr['ext_start'] = onehit.best_start
        rr['ext_end'] = onehit.best_end
        rr['blast_text'] = blasthsp2pre(onehit.source.annotations['blast'][1])
        rr['pictures'] = _prepare_pictures(ext)
        rr['eval'] = onehit.source.annotations['blast'][1].expect

        # create seqviewurl here
        es = onehit.source.annotations['extended_start']
        ee = onehit.source.annotations['extended_end']
        if es > ee:
            es, ee = [ee, es]

        diff = 1000 + 2*len(ext)
        es -= diff
        ee += diff

        if getattr(data.args, 'show_gene_browser', False):
            rr['draw_seqview'] = True
        else:
            rr['draw_seqview'] = False

        rr['seqviewurl'] = ''.join([
            '?embeded=true&',
            'noviewheader=true&',
            'id={}&'.format(onehit.source.annotations['blast'][0]),
            'v={}:{}&'.format(es, ee),
            'appname=rna_blast_analyze_beta&',
            'mk={}:{}|BestMatch!&'.format(onehit.best_start, onehit.best_end),
            'multipanel=false',
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
            # a = f.read()
            # np['pic'] = "data:image/svg;base64," + base64.b64encode(f.read().encode()).decode()
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
        command = ''

    # parameters
    params = [i for i in dir(data.args) if not i.startswith('__') and not callable(getattr(data.args, i))]
    p_text = []
    for arg in params:
        p_text.append((arg, getattr(data.args, arg)))

    if not hasattr(data, 'date_of_run'):
        data.date_of_run = None

    prepared_footer_data = {
        'command': command,
        'parameters': p_text,
        'exec_date': data.date_of_run,
    }
    return prepared_footer_data
