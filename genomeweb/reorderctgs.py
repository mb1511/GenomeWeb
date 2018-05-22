'''
Created on 15 May 2018

@author: mb1511
@organization: University of Bristol
@contact: mb1511@bristol.ac.uk
@summary:
'''
from __future__ import print_function
from builtins import range
import re

from genomeweb import fasta
from genomeweb.blast import local_blast
from os.path import join

def _get_props(gene_name):
    p = dict([ x.split('=') for x in re.findall('(?<=\[)(.*?)(?=\])', gene_name) if len(x.split('=')) == 2])
    contig = re.findall('(?<=\|)(.*?)(?=_cds)', gene_name)
    if not contig:
        contig = re.findall('(?<=\>)(.*?)(?=_cds)', gene_name)
    try:
        p['contig'] = contig[0]
    except IndexError:
        # no contig
        pass 
    return p

def _get_loc(props):
    if 'comp' in props['location']:
        props['location'] = props['location'][11:-1]
    if 'join' in props['location']:
        props['location'] = props['location'][5:-1].split(',')
        props['location'] = props['location'][0] + '..' + props['location'][1]
    l =  props['location'].split('..')
    f = re.sub('[<;>]|&lt|&gt', '', l[0])
    t = re.sub('[<;>]|&lt|&gt', '', l[1])
    return (f, t)

def run(
        query, reference, chunk=1000, step=2000, order_out='order_out.fna',
        wd='', run_short=True, shortck=30, shortsp=100, l_max_hsps=4, **kw):
    '''
    TODO: add doc string
    '''
    
    # split genome into bits
    n1 = fasta.fasta_read(query, generator=False)    # query
    n2 = fasta.fasta_read(reference, generator=False)    # reference
    
    ws_1, ws_2 = [], []
    nams1, nams2 = [], []
    full = False
    short = False
    
    with open(join(wd, 'nuc_1.fna'), 'w') as s, \
         open(join(wd, 'nuc_short.fna'), 'w') as d:
        
        for num, ctg in enumerate(n1):
            if len(ctg) < chunk and run_short:
                # if short contig, run short analysis
                short = True
                ck = shortck
                sp = shortsp
            else:
                full = True
                ck = chunk
                sp = step
            
            for i in range(0, len(ctg), sp):  #chunk
                if i + ck <= len(ctg):
                    seq = ctg[i : i + ck]
                    name = re.sub('[\n\r]', '', ctg.name)
                    
                    name += '_[ctg=%d] [location=%d..%d]\n' % (
                        num,
                        i+1,
                        1+i+len(seq))
                    
                    if ck == chunk:
                        s.write(name)
                        s.write(seq + '\n')
                        full += 1
                    else:
                        d.write(name)
                        d.write(seq + '\n')
                        short = True
            
            ws_1.append(len(ctg))
            nams1.append(ctg.name[1:ctg.name.find(' ')])
    
    with open(join(wd, 'nuc_2.fna'), 'w') as s:
        for i, ctg in enumerate(n2):
            name = re.sub('[\n\r]', '', ctg.name) + '_[ctg=%d]\n' % i
            s.write(name)
            s.write(ctg.seq)
            s.write('\n')
            ws_2.append(len(ctg))
            nams2.append(ctg.name[1:ctg.name.find(' ')])
    
    # set blast defaults - can be overwritten in **kw
    short_defaults = dict(
        db_path=join(wd, 'db.fna'),
        mr=0, mev=10000,
        join_hsps=False,
        b_type='blastn',
        r_a=True,
        num_threads=4)
    short_defaults.update(kw)
    
    if 'max_hsps' in kw:
        del kw['max_hsps']
    
    long_defaults = dict(make_db=False)
    
    defaults = dict(
        db_path=join(wd, 'db.fna'),
        mr=0, mev=10000,
        join_hsps=False,
        b_type='blastn',
        r_a=True,
        num_threads=4,
        max_hsps=l_max_hsps)
    
    defaults.update(kw)
    long_defaults.update(defaults)
    
    if short:
    
        local_blast.run(
            join(wd, 'nuc_2.fna'),
            query=join(wd, 'nuc_short.fna'),
            out=join(wd, 'temp_blast_1.xml'),
            outfmt='details', task='blastn-short',
            **short_defaults)
        if full:
            local_blast.run(
                join(wd, 'nuc_2.fna'),
                query=join(wd, 'nuc_1.fna'),
                out=join(wd, 'temp_blast_2.xml'),
                outfmt='details',
                **long_defaults)
        else:
            with open(join(wd, 'temp_blast_2.xml'), 'w') as fbx:
                # clear file contents
                fbx.write('')
            
        # concatenate blast output files
        with open(join(wd, 'temp_blast.xml'), 'w') as o, \
             open(join(wd, 'temp_blast_1.xml')) as i1,   \
             open(join(wd, 'temp_blast_2.xml')) as i2:
            
            for line in i1:
                if '</BlastOutput_iterations>' not in line:
                    o.write(line)
                else:
                    break
            
            read_line=False
            for line in i2:
                if '<Iteration>' in line or read_line:
                    read_line = True
                    o.write(line)
            # read concatenated xml output but don't run
            itrs = local_blast.run(
                join(wd, 'nuc_2.fna'),
                query=join(wd, 'nuc_1.fna'),
                out=join(wd, 'temp_blast.xml'),
                mr=0, mev=10000,
                outfmt='details', join_hsps=False,
                b_type='blastn', r_a=True,
                blast_run=False, make_db=False)
    
    else:
        itrs = local_blast.run(
            db=join(wd, 'nuc_2.fna'),
            query=join(wd, 'nuc_1.fna'),
            out=join(wd, 'temp_blast.xml'),
            outfmt='details',
            **defaults)
    
    
    def _avg(lst):
        out = []
        tot = len(lst)
        n = len(lst[0])
        bin_width = 50000
        bin_bounds = range(0, 10000000, bin_width)
        bins = [0 for _ in bin_bounds]
        for i in range(n):
            if i:
                s = sum(float(x[i]) for x in lst)
                out.append(int(s/tot))
            else:
                regions = [x[i] for x in lst]
                for r in regions:
                    # find and add 1 to bin
                    for bi, b in enumerate(bin_bounds):
                        if r >= b and r < b + bin_width:
                            bins[bi] += 1
                            break
                out.append(bin_bounds[max(range(len(bins)), key=lambda x: bins[x])])               
        return tuple(out)
    
    ctg_pos = {}
    
    for alns in itrs:
        for hit in alns:
            
            # get ctg index and posistion
            q_props = _get_props(hit.q_def)           
            q_ctg_num = int(q_props['ctg'])
            #q_loc = int(_get_loc(q_props)[0])
                 
            # get reference ctg index and position
            h_props = _get_props(hit.h_def)
            r_ctg_num = int(h_props['ctg']) # should be 0 for reference of length 1
            hit_pos_array = []
            if len(hit.h_from) > 0:
                for i, hsp_loc in enumerate(hit.h_from):
                    h_loc = int(hsp_loc)
                    h_to = int(hit.h_to[i])
                    h_len = abs(h_loc - h_to)
                    if h_loc < h_to:
                        direction = 1
                    else:
                        direction = -1
                    pid = hit.pid[i]
                    hit_pos_array.append((h_loc, r_ctg_num, h_len, pid, direction))
            else: 
                continue
            hit_pos_array = sorted(hit_pos_array, key=lambda x: (-x[1], x[2], x[3]), reverse=True)
            try:
                ctg_pos[q_ctg_num]
                ctg_pos[q_ctg_num].append(hit_pos_array[0])
            except KeyError:
                ctg_pos[q_ctg_num] = [hit_pos_array[0],]
    order_ctgs = []
    for k, i in ctg_pos.items():
        # take average of (up to) top 50 positions
        order_ctgs.append((k, _avg(sorted(i, key=lambda x: (-x[1], x[2], x[3]), reverse=True)[:50])))
        #print(order_ctgs[-1])
    order_ctgs = sorted(order_ctgs, key=lambda x: (x[1][1], x[1][0]))
    with open(order_out, 'w') as w:
        for i in order_ctgs:
            #print(i)
            w.write(n1[i[0]].name + '\n')
            if i[1][4] < 0:
                w.write(n1[i[0]].reverse_complement + '\n')
            else:
                w.write(n1[i[0]].seq + '\n')
    if 'quiet' in kw:
        quiet = kw['quiet']
    else:
        quiet = False
    if not quiet:   
        print('Contig reorder complete. Output file: %s.' % order_out)
    return(order_out)
    