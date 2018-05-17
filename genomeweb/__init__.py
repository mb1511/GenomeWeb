'''
Created on 6 Feb 2017

@author: mb1511
'''

from __future__ import print_function

from pyveplot2 import Hiveplot, Axis, Node
import svgwrite as sw
import re
from os.path import basename, splitext, join
import glob
import numpy as np
import colorcet as cc

import fasta
from blast import local_blast
import reorderctgs

def _get_matches(
        n1, n2, chunk=2000, step=5000, pid_cov=90, hit_len=1500,
        short_ctgs=False, shortck=30, shortsp=100, wd=''):
    
    # should return ctg num and match posistion for n1 and n2
    # split genome into bits
    n1 = fasta.fasta_read(n1, generator=False)    # genome 1
    n2 = fasta.fasta_read(n2, generator=False)    # genome 2
    
    ws_1, ws_2 = [], []
    nams1, nams2 = [], []
    full=0
    short = False
    
    with open(join(wd, 'nuc_1.fna'), 'w') as s, \
         open(join(wd, 'nuc_short.fna'), 'w') as d:
        
        for num, ctg in enumerate(n1):
            if short_ctgs and len(ctg) < chunk:
                # if short contig, run short analysis
                ck = shortck
                sp = shortsp
            else:
                ck = chunk
                sp = step
            
            for i in xrange(0, len(ctg), sp):  #chunk
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
    
    if short:
    
        local_blast.run(
            join(wd, 'nuc_2.fna'),
            db_path=join(wd, 'db.fna'),
            query=join(wd, 'nuc_short.fna'),
            out=join(wd, 'temp_blast_1.xml'),
            mr=0, mev=10000,
            outfmt='details', join_hsps=False,
            b_type='blastn', r_a=True,
            num_threads=4,task='blastn-short')
        if full:
            local_blast.run(
                join(wd, 'nuc_2.fna'),
                db_path=join(wd, 'db.fna'),
                query=join(wd, 'nuc_1.fna'),
                out=join(wd, 'temp_blast_2.xml'),
                mr=0, mev=10000,
                outfmt='details', join_hsps=False,
                b_type='blastn', r_a=True,
                num_threads=4, max_hsps=4, make_db=False)
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
                db_path=join(wd, 'db.fna'),
                query=join(wd, 'nuc_1.fna'),
                out=join(wd, 'temp_blast.xml'),
                mr=0, mev=10000,
                outfmt='details', join_hsps=False,
                b_type='blastn', r_a=True,
                blast_run=False, make_db=False)
    
    else:
        itrs = local_blast.run(
            db=join(wd, 'nuc_2.fna'),
            db_path=join(wd, 'db.fna'),
            query=join(wd, 'nuc_1.fna'),
            out=join(wd, 'temp_blast.xml'),
            mr=0, mev=10000,
            outfmt='details', join_hsps=False,
            db_type='nucl',
            b_type='blastn', r_a=True,
            num_threads=4, max_hsps=4)
    
    print('BLAST Search Complete.')    
    
    for alns in itrs:
        for hit in alns:
            
            q_props = _get_props(hit.q_def)
            
            # get axis index
            l_axis_index = int(q_props['ctg'])
            # get location on (left) contig
            q_loc = int(_get_loc(q_props)[0])
            
            # get right contig y coordinates
            h_props = _get_props(hit.h_def)
            r_axis_index = int(h_props['ctg'])
            
            if len(hit.h_from) > 0:
                for i, hsp_loc in enumerate(hit.h_from):
                    h_loc = int(hsp_loc)
                    if q_props['ctg'] == '-1':
                        pid = 100 * float(hit.idnt[i])/float(hit.q_len)
                    else:
                        pid = hit.pid[i]
                    if pid >= pid_cov and abs(int(hit.h_to[i]) - h_loc) >= hit_len:    
                        yield (q_loc, l_axis_index, h_loc, r_axis_index, pid)
    
def create_web(
        genome_array=[], reference_genome='',
        working_directory='', out_file='default.svg',
        include_reference=False,
        reorder=True,
        palette='bgy',
        palette_usage=0.8,
        bezier_max_n=4
        connection_opts=dict(
            stroke_width='0.34', stroke_opacity='0.4'),
        axes_opts=dict(
            stroke='black', stroke_width='1'),
        inner_radius=100,
        outer_radius=300,
        boarder_offset=80,
        x_scaling='200%',
        y_scaling='200%',
        add_labels=True,
        label_offset=1.1,
        reorder_opts=dict(),
        matches_opts=dict()):
    '''
    TODO: add doc string
    '''
    
    if include_reference:
        genome_array.append(reference_genome)
    
    print('Number of Genomes: ', len(genome_array))
    
    names = [splitext(basename(g))[0] for g in genome_array]
    
    # reindex all contigs against reference
    # -> new adjusted array
    if reorder:
        for i, genome in enumerate(genome_array):
            if not include_reference or i != len(genome_array)-1:
                genome_array[i] = reorderctgs.run(
                    genome, reference_genome,
                    order_out=join(
                        working_directory,
                        '%s_reorder.fna' % splitext(basename(genome))[0]),
                    make_db=not i,
                    wd=working_directory,
                    **reorder_opts)
    
    matches = []
    contig_sizes = []
    genome_sizes = []
    # get matching regions for adjacent pairs in list starting with [-1] v [0]
    # line filtering is alos done here
    for i, g in enumerate(genome_array):
        # get genome contig lengths
        ctgs = fasta.fasta_read(g)
        contig_sizes.append([len(ctg) for ctg in ctgs])
        genome_sizes.append(sum(contig_sizes[-1]))     
        # adds generator to matches list
        matches.append(_get_matches(
            genome_array[i-1], g, wd=working_directory **matches_opts))
        # will return ctg_l, pos_l, ctg_r, pos_r and pid
        # when evaluated
    
    # ======== Define Axes ========
    
    # inner-shell radius
    pre_r = inner_radius
    # outer-shell radius
    r = outer_radius + pre_r
    # outer-shell boarder offset
    offset = boarder_offset
    # origin    
    ox = r + pre_r + offset
    oy = r + pre_r + offset
    positions = []
    theta = (2 * np.pi)/len(genome_array)
    for i in xrange(len(genome_array)):
        positions.append((
            (int(ox + pre_r * np.cos(theta * i)), int(oy + pre_r * np.sin(theta * i))),
            (int(ox + r * np.cos(theta * i)), int(oy + r * np.sin(theta * i)))
            ))
    hive = Hiveplot(out_file, size=(x_scaling, y_scaling))
    hive.axes = [Axis(
        start=a, end=b,
        **axes_opts
        ) for a, b in positions]
    
    # ======== Connect up the Axes ========
    
    if 'pid_cov' in matches_opts:
        pid_cov = matches_opts['pid_cov']
    else:
        pid_cov = 90.0
    n_paths = 0
    for i, match in enumerate(matches):
        for node_index, (pos1, ctg1, pos2, ctg2, pid) in enumerate(match):
            sum1 = pos1
            for j, s in enumerate(contig_sizes[i-1]):
                if j < ctg1:
                    sum1 += s
            off1 = float(sum1)/float(genome_sizes[i-1])
            hive.axes[i-1].add_node(Node(), offset=off1)
            
            sum2 = pos2
            for j, s in enumerate(contig_sizes[i]):
                if j < ctg2:
                    sum2 += s
            off2 = float(sum2)/float(genome_sizes[i])
            hive.axes[i].add_node(Node(), offset=off2)
            
            # change this depending on color scheme used
            if pid < pid_cov:
                pid = pid_cov
            pid = int(255 * (pid - pid_cov) * palette_usage)
            n_paths += 1
            hive.connect(
                hive.axes[i-1], node_index, np.pi/4,
                hive.axes[i], node_index, -np.pi/4,
                curved=False if len(genome_array) > bezier_max_n else True,
                stroke=cc.palette[palette][pid],
                **connection_opts)
        
        # Clear node arrays
        hive.axes[i-1].nodes = []
        hive.axes[i].nodes = []
    
    # ======== Add Labels to Axes ========
    
    if add_labels:
        for i, g in enumerate(names):
            hive.axes[i].add_node(Node(g, draw_label=True), offset=label_offset)
    
    print('Total number of connections: ', n_paths)
    # Save and write svg
    print('Saving %s...' % out_file)
    hive.save()
    print('Finished.')
    return 0