'''
Created on 6 Feb 2017

@author: mb1511
'''

from __future__ import print_function

from pyveplot2 import Hiveplot, Axis, Node
import svgwrite as sw
import re
from os.path import basename, splitext, join, exists
import numpy as np
import colorcet as cc

import fasta
from blast import local_blast
import reorderctgs

def _get_props(gene_name):
    p = dict([ x.split('=') for x in re.findall('(?<=\[)(.*?)(?=\])', gene_name) if len(x.split('=')) == 2 ])
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
        palette_usage=1.0,
        bezier_max_n=4,
        connection_opts=dict(
            stroke_width='0.34', stroke_opacity='0.4'),
        axes_opts=dict(
            stroke='black', stroke_width='1'),
        inner_radius=30,
        outer_radius=140,
        border_offset=10,
        x_scaling='500px',
        y_scaling='500px',
        add_labels=True,
        label_offset=1.1,
        reorder_opts=dict(),
        matches_opts=dict(),
        svg_opts=dict()):
    '''
    Create Gemoic Comparison Web
    
    Usage:
    
    >>> genomeweb.create_web(genome_list, reference_genome, **options)
    
    Necessary Arguments:
        
        Name               Type     Description
        
        genome_array       list     list of paths to genomes                
        reference_genome   str      path to reference to order contigs
                                    by (not required if reorder=False)

    Standard Options:

        working_directory  str      path to scratch space for program
                                    to write files
        out_file           str      path to output SVG file
        include_referecne  bool     include referecne genome in
                                    resulting map                            
        reorder            bool     reorder contigs against reference
        palette            str      color palette to use from 
                                    colorcet.palette
        add_labels         bool     add labels to axes
        label_offset       float    label offset
        inner_radius       int      inner-shell radius
        outer_radius       int      distance from inner shell to outer-
                                    shell bounds
        border_offset      int      distance from border (increase to
                                    correctly display longer genome
                                    names)
        matches_opts       dict     kwargs to parse to find matches 
                                    function, such as %id cutoff value
                                    (pid_cov) and minimum hit length
                                    (hit_len) for use in filtering hits
                                
    Advanced Options:

        palette_usage      float    decimal percent of palette spectrum
                                    to use
        bezier_max_n       int      max number of genomes before
                                    straight lines are used instead of
                                    bezier curves (set to 0 for always
                                    straight)
        x_scaling          str      SVG x scale factor
        y_scaling          str      SVG y scale factor
        connection_opts    dict     dictionary of line options
        axes_opts          dict     dictionary of axes options
        reorder_opts       dict     kwargs to parse to reorder function
        svg_opts           dict     additional properties for base SVG
                                    (see svgwrite for docs)
      
    '''

    # check all files exist before running anything
    if isinstance(reference_genome, (list, tuple)):
        # check the length of the two lists match if using multiple
        # reference genomes
        assert len(reference_genome) == len(genome_array), 'Unequal number of references and query genomes.'
        for f in reference_genome:
            if f:
                assert exists(f), '%s cannot be found.' % f
    else:
        assert exists(reference_genome), '%s cannot be found.' % reference_genome
    for f in genome_array:
        assert exists(f), '%s cannot be found.' % f
        
    if include_reference:
        genome_array.append(reference_genome)
    
    print('Number of Genomes: ', len(genome_array))
    
    if isinstance(reference_genome, (list, tuple)):
        print('Number of unique reference genomes: %d' % len(set(reference_genome)))
    
    names = [splitext(basename(g))[0] for g in genome_array]
    
    # reindex all contigs against reference
    # -> new adjusted array
    if reorder:
        ref_to_use = ''
        for i, genome in enumerate(genome_array):
            if not include_reference or i != len(genome_array)-1:
                if isinstance(reference_genome, (list, tuple)):
                    ref = reference_genome[i]
                    p_ref = reference_genome[-1]
                    if ref:
                        ref_to_use = ref
                    # only make new db for unique/non-empty/initial
                    # genome in referecne list, else use previous entry
                    genome_array[i] = reorderctgs.run(
                        genome, ref_to_use,
                        order_out=join(
                            working_directory,
                            '%s_reorder.fna' % splitext(basename(genome))[0]),
                        make_db=(not i or ref != p_ref) and ref,
                        wd=working_directory,
                        **reorder_opts)
                else:
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
            genome_array[i-1], g, wd=working_directory, **matches_opts))
        # will return ctg_l, pos_l, ctg_r, pos_r and pid
        # when evaluated
    
    # ======== Define Axes ========
    
    # inner-shell radius
    pre_r = inner_radius
    # outer-shell radius
    r = outer_radius + pre_r
    # outer-shell boarder offset
    offset = border_offset
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
    xb = int(re.search(r'\d+', x_scaling).group())
    yb = int(re.search(r'\d+', y_scaling).group())
    hive = Hiveplot(
        out_file,
        size=(x_scaling, y_scaling),
        viewBox='0 0 %d %d' % (xb, yb),
        **svg_opts)
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
            pid = int(255 * palette_usage * (pid - pid_cov) / (100 - pid_cov) )
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