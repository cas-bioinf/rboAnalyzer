import argparse


def db2shape(seq):
    """converts Vienna dot bracket representation of rna structure to a shape level 5"""
    [c, nest] = nesting(seq, ('(', ')'))
    [b, s, rec_level] = analyze(nest, 0)
    if b == rec_level:
        # means that no basepairs where found and b, rec_level are empty lists
        b = ['_']
    return b, s, rec_level


def nesting(seq, br=('(',')'), gap_mark=0, initial_val=1, input_gap_mark='.'):
    """br must be nesting brackets of the kind you want to analyze"""
    c = initial_val
    nest = []   # nest should be same length as seq at the end
    for char in seq:
        if char == br[0]:
            nest.append(c)
            c += 1
        elif char == br[1]:
            c -= 1
            nest.append(c)
        elif char == input_gap_mark:
            nest.append(gap_mark)
        else:
            raise ValueError('Encountered char, which is not present in input character set: char:{}'
                             ' input brackets: {}{}, input gap character: {}'.format(char,
                                                                                     br[0],
                                                                                     br[1],
                                                                                     input_gap_mark))
    return c, nest


def analyze(n, rec):
    rec += 1
    c = -1
    #   print('rec entry {}'.format(rec))
    #   print(n)
    lb = []
    rb = []
    qq = []
    rec_level = []
    # copy n
    nn = n.copy()
    nn.reverse()
    maxn = max(n)
    bracket = False
    while True:
        c += 1
    # for c in range(len(n)):
        if c >= len(n):
            break
        if n[c] <= 0:
            continue

        if n[c] > 0 and n.count(n[c]) > 2:
            #   print('branch at {} c pos {}'.format(n[c], c))
            if bracket:
                # add bracket if bracket was present in stack
                lb.append('[')
                rb.append(']')
            # todo add from which point to which point it should analyze and return ending point of analysis and shift c
            # todo this will be executed always, if there is any bracket
            start = n.index(n[c])
            end = len(nn) - nn.index(n[c])
            splited = _split_stems(n[start:end])
            for q in splited:
                [sub_sh, r, rr] = analyze(q, rec)
                qq += sub_sh
                c += len(q)
                rec_level += rr

        else:
            bracket = True
            if n[c] == maxn:
                # midle of the stack
                # if rec == 1:
                #     pass
                # else:
                lb.append('[')
                rb.append(']')
                break

    # consolidate
    sh = lb + qq + rb
    if bracket:
        rec_level = _prepend(rec_level, rec)
        rec_level.append(rec)
    return sh, rec, rec_level


def _prepend(s,el):
    # prepends element to a list
    s.reverse()
    s.append(el)
    s.reverse()
    return s


def _split_stems(nest):
    """split stems denoted in nesting vector and return spliced vectors with indices by one lower"""
    # stems
    st = [i for i, x in enumerate(nest) if x == nest[0]]
    st.append(len(nest))
    qq = []
    for i in range(0,len(st)-1,2):
        r = nest[st[i]:st[i+2]]
        qq.append(r)
    return qq


def f_arg():
    parser = argparse.ArgumentParser(
        description='Dot bracket to shape notation function\n'
                    'In order to replace slow RNAshapes abstract\n'
                    'do not work on pseudoknots !\n'
                    'if input structure is broken i.e. not correctly bracketed,'
                    ' it will end in error\n'
                    'execution: Function will return [] for continuous basepair stem,'
                    ' ignoring bulges. The unpaired sections are not reported.'
                    ' Branching is recognized in the same manner i.e. ignoring bulges.\n'
                    'For demonstration run db2shape -demo'
    )
    parser.add_argument(
        '-demo', action='store_true', help='run demo'
    )
    parser.add_argument(
        '-bench', default=None, help='compare with Rapidshapes level 5, input is file like for rapidshapes output is to inputfile+_db2shape'
    )
    parser.add_argument(
        '-str', default=None, help='structure in Vienna format to convert to shape level 5,'
                                   ' lonely pairs are not allowed '
                                   'for the sake of reproducibility of RNAshapes results'
    )
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    import os
    from tempfile import mkstemp
    import re

    args = f_arg()
    if args.str:
        [a,b,c] = db2shape(args.str)
        print('>shape level 5\n')
        print(*a)
        print(*c)

    if args.demo:
        l = []
        l.append(
            '..(((((((((.((((((((((((((((((((((....))))))))))....'
            '........))))))))).))).............(((((((....)))))))'
            '((((((((((...))))))))))....................)))))))))'
            '...........(((((((((((....)))))))))))......'
        )
        # sh1: _[_[[[]_]_]_[][]_]_[]_
        # sh2: [[[[]_]_][][]][]
        # sh3: [[[[]]][][]][]
        # sh4: [[][][]][]
        # sh5: [[][][]][]

        l.append(
            '(((((((((....((((((((..(((..((((..((....))..))))..)))'
            '...)))))).))(((((((.....(((((.(((....))))))))....))))'
            ')))))))))))).'
        )
        # 5 [[][]]
        # 4 [[[[[]]]][[]]]
        # 3 [[[[[[]]]]][[[]]]]
        # 2 [[[_[_[_[]_]_]_]_][_[_[]]_]]
        # 1 [_[[_[_[_[]_]_]_]_][_[_[]]_]]_

        l.append('(((())(())))')
        for s in l:
            [a,b,c] = db2shape(s)
            print('Structure:')
            print(s)
            print('Shape')
            print(*a)
            print('nesting level')
            print(*[str(i) for i in c])

    if args.bench:
        # compute db2shape
        # read the file first is header, second is structure sequence in Vienna dot bracket format
        fp, filename = mkstemp(prefix='rba_', suffix='_16', dir=CONFIG.tmpdir)
        ss = []
        with open(args.bench, 'r') as f, os.fdopen(fp, 'w') as fw:
            txt = f.readline()
            while txt:
                if txt.strip() == '':
                    continue
                txt = txt.strip()
                if txt[0] == '>':
                    # parse
                    txt2 = f.readline()
                    if '(' and ')' and '.' not in txt2:
                        txt2 = f.readline()
                    structure = re.search('[\(\.\)]+', txt2).group()
                    ss.append([txt, structure])
                    # write for safety
                    fw.write(txt + '\n' + structure + '\n')
                txt = f.readline()

        # actualy compute dbshape
        answer = []
        for s in ss:
            [a, b, c] = db2shape(s[1])
            answer.append(''.join(a))

        # insert here ane output file
        with open(args.bench + '_db2shape_test.txt', 'w') as fout:
            for a in answer:
                fout.write(a + '\n')
