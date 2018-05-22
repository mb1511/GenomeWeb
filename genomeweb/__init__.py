'''
Created on 6 Feb 2017

@author: mb1511
'''

from __future__ import print_function
from builtins import range
from genomeweb.pyveplot2 import Hiveplot, Axis, Node
import svgwrite as sw
import svgutils as su
import re
from os.path import basename, splitext, join, exists
import numpy as np
import colorcet as cc

from genomeweb import fasta, reorderctgs
from genomeweb.blast import local_blast
#import reorderctgs

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
        short_ctgs=False, shortck=30, shortsp=100, wd='', quiet=False):
    
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
    
    if short:
    
        local_blast.run(
            join(wd, 'nuc_2.fna'),
            db_path=join(wd, 'db.fna'),
            query=join(wd, 'nuc_short.fna'),
            out=join(wd, 'temp_blast_1.xml'),
            mr=0, mev=10000,
            outfmt='details', join_hsps=False,
            b_type='blastn', r_a=True,
            num_threads=4,task='blastn-short',
            quiet=quiet)
        if full:
            local_blast.run(
                join(wd, 'nuc_2.fna'),
                db_path=join(wd, 'db.fna'),
                query=join(wd, 'nuc_1.fna'),
                out=join(wd, 'temp_blast_2.xml'),
                mr=0, mev=10000,
                outfmt='details', join_hsps=False,
                b_type='blastn', r_a=True,
                num_threads=4, max_hsps=4, make_db=False,
                quiet=quiet)
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
                blast_run=False, make_db=False,
                quiet=quiet)
    
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
            num_threads=4, max_hsps=4,
            quiet=quiet)
    if not quiet:
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
        include_reference=False, reorder=True,
        palette='bgy', palette_usage=1.0,
        bezier_max_n=4,
        connection_opts=dict(
            stroke_width='0.34', stroke_opacity='0.4'),
        axes_opts=dict(
            stroke='black', stroke_width='1'),
        size=500, inner_radius=12, outer_radius=64, border_offset=None,
        rotation=0, x=0, y=0, width=None, height=None, viewBox=None,
        label_names=[], add_labels=True, label_offset=1.1,
        font_size=11, font_family='Arial', custom_font=None,
        reorder_opts=dict(),
        matches_opts=dict(),
        svg_opts=dict(),
        append=False, quiet=False):
    '''
    Create Gemoic Comparison Web
    
    Necessary Arguments:
        
        Name               Type     Description
        
        genome_array       list     list of paths to genomes                

    Standard Options:
    
        reference_genome   str      path to reference to order contigs
                                    by (required if reorder=True)
        working_directory  str      path to scratch space for program
                                    to write files
        out_file           str      path to output SVG file
        include_referecne  bool     include referecne genome in
                                    resulting map                            
        reorder            bool     reorder contigs against reference
        palette            str      color palette to use from 
                                    colorcet.palette
        add_labels         bool     add labels to axes
        label_offset       float    label offset, 1.0-1.5 should be
                                    fine, default = 1.1
        inner_radius       int      inner-shell radius (% of size) 
        outer_radius       int      inner shell to outer-shell bounds
                                    (% of size)
        matches_opts       dict     kwargs to parse to find matches 
                                    function, such as %id cutoff value
                                    (pid_cov) and minimum hit length
                                    (hit_len) for use in filtering hits
                                    
    Basic geometry options
    
        size               int      size of hive plot in px
        x                  int      x origin, default = 0
        y                  int      y origin, dafault = 0
        rotation           float    value in degrees to rotate all axes
                                    by (clockwise direction)
                                    
    Mulit-panel geometry options:
        
        append             bool     append output to existing SVG,
                                    default to append is out_file
        append             str      file to append
        width              int      width of final SVG in px
        height             int      height of fnial SVG in px
                        
    Advanced Options:

        palette_usage      float    decimal percent of palette spectrum
                                    to use
        bezier_max_n       int      max number of genomes before
                                    straight lines are used instead of
                                    bezier curves (set to 0 for always
                                    straight)
        custom_font        str      full custom SVG text element for
                                    labels e.g.
                                    "font-size:12px; font-family:Arial;
                                    font-style:italic"
        connection_opts    dict     dictionary of line options
        axes_opts          dict     dictionary of axes options
        reorder_opts       dict     kwargs to parse to reorder function
        svg_opts           dict     additional properties for base SVG
                                    (see svgwrite for docs)
    
    Additional options (likely do not need to be altered):
    
        border_offset      int      distance from border, default
                                    calculated automatically
        viewBox            tuple    viewBox for SVG, default
                                    calculated automatically        
    
    Usage:
    
        >>> fig = genomeweb.create_web(genome_list, **options)
        >>> # Returns svgutils.transform.SVGFigure for any extra post editing.
        >>> # e.g. add an "A" to the top left of the figure.
        >>> import svgutils.transform as sg
        >>> fig.append(sg.TextElement(0, 20, 'A', size=12))
        >>> fig.save('figure.svg')
    
    '''

    # check all files exist before running anything
    if isinstance(reference_genome, (list, tuple)):
        # check the length of the two lists match if using multiple
        # reference genomes
        assert len(reference_genome) == len(genome_array), 'Unequal \
number of references and query genomes.'
        for f in reference_genome:
            if f:
                assert exists(f), '%s cannot be found.' % f
    elif reorder:
        assert exists(reference_genome), 'Reference genome: %s cannot be found.' % reference_genome
    for f in genome_array:
        assert exists(f), '%s cannot be found.' % f
        
    if include_reference:
        genome_array.append(reference_genome)
    if not quiet:
        print('Number of Genomes: ', len(genome_array))
    
    if isinstance(reference_genome, (list, tuple)):
        if '' not in reference_genome:
            unq = len(set(reference_genome))
        else:
            unq = len(set(reference_genome)) - 1
        if not quiet:
            print('Number of unique reference genomes: %d' % unq)
            
    
    # ======== Setup and Check Geometry ========
    
    # inner-shell radius
    pre_r = int(inner_radius * size / 200)
    # outer-shell radius
    r = int(outer_radius * size / 200) + pre_r
    if border_offset:
        # outer-shell boarder offset
        offset = int(border_offset * size / 200)
    else:
        offset = size / 2 - r
        assert offset >= 0, 'Some axes may not be visible in resulting\
 output. Decrease radii percentages to below 100 total.'
    assert r + offset <= size / 2, 'Some axes may not be visible in \
resulting output. Decrease radii percentages to below 100 total.'
    # origin    
    ox = x + r + offset
    oy = y + r + offset
    
    # ======== Re-index Contigs ========
    
    # set label names
    if not label_names:
        names = [splitext(basename(g))[0] for g in genome_array]
    else:
        assert len(label_names) >= len(genome_array), 'Number of labels\
 does not match number of genomes.'
        names = label_names
    
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
    
    # axis geometry
    positions = []
    if genome_array:
        n_axes = len(genome_array)
    else:
        n_axes = len(label_names)
    theta = (2 * np.pi)/n_axes
    i_theta = np.deg2rad(rotation)
    for i in range(n_axes):
        angle = i_theta + theta * i
        positions.append((
            (
                int(ox + pre_r * np.cos(angle)),
                int(oy + pre_r * np.sin(angle))),
            (
                int(ox + r * np.cos(angle)),
                int(oy + r * np.sin(angle)))
            ))
    
    # set global geometry
    if width is None:
        width = size
    if height is None:
        height = size
    if viewBox is None:
        viewBox = (0, 0, width, height)
        
    # save previous SVG in memory, if in append mode
    if append:
        if isinstance(append, bool):
            template = su.transform.fromfile(out_file)
        else:
            template = su.transform.fromfile(append)
    
    # define hive plot object
    hive = Hiveplot(
        out_file,
        size=('%dpx' % width, '%dpx' % height),
        viewBox='%d %d %d %d' % viewBox,
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
    if not custom_font:
        font='font-size:%dpx; font-family:%s' % (font_size, font_family)
    else:
        font=custom_font
    if add_labels:
        for i, g in enumerate(names):
            hive.axes[i].add_node(
                Node(g, draw_label=True),
                offset=label_offset, font=font)
    
    # ======== Save and write svg ========
    if not quiet:
        print('Total number of connections: ', n_paths)
        if append:
            print('Saving and merging files...')
        else:
            print('Saving %s...' % out_file)
    hive.save()
    if append:
        template.append(su.transform.fromfile(out_file))
        template.save(out_file)
    else:
        template = su.transform.fromfile(out_file)
    if not quiet:
        print('Finished.')
    return template