'''
Created on 6 Feb 2017

@author: mb1511
'''
from __future__ import print_function
from builtins import range

import re
import numpy as np
import logging
from warnings import warn
from os.path import basename, splitext, join, exists

import svgwrite as sw
import svgutils as su
import colorcet as cc

from genomeweb import fasta, reorderctgs
from genomeweb.pyveplot2 import Hiveplot, Axis, Node
from genomeweb.blast import local_blast, match_calc
from genomeweb.reorderctgs import _get_props, _get_loc

class GeometryWarning(UserWarning):
    pass

def _get_matches(
        n1_path, n2_path, chunk=2000, step=5000, pid_cov=90, hit_len=1500,
        short_ctgs=False, shortck=30, shortsp=100, wd='',
        l_max_hsps=4, quiet=False, **kw):
    
    # should return ctg num and match posistion for n1 and n2
    # split genome into bits
    n1 = fasta.fasta_read(n1_path, generator=False)
    n2 = fasta.fasta_read(n2_path, generator=False)
    
    itrs = match_calc.get_matches(
        n1, n2, n1_path, n2_path, chunk, step, short_ctgs, shortck, shortsp,
        wd, l_max_hsps, quiet, **kw)  
    
    if not quiet:
        print('BLAST Search Complete.')    
    
    for alns in itrs:
        for hit in alns:
            
            q_props = _get_props(hit.q_def)
            
            # get axis index
            l_axis_index = int(q_props['ctg'])
            # get location on contig
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
        palette='bgy', palette_usage=1.0, invert=False, flip=False,
        connection_opts=dict(
            stroke_width='0.34', stroke_opacity='0.4'),
        axes_opts=dict(
            stroke='black', stroke_width='1'),
        size=500, inner_radius=12, outer_radius=64, border_offset=None,
        rotation=0, x=0, y=0, width=None, height=None, viewBox=None,
        label_names=[], add_labels=True, label_offset=1.1,
        font_size=11, font_family='Arial', custom_font=None,
        bezier_max_n=4, source_angle=None, target_angle=None, dummy_axes=0,
        reorder_opts=dict(),
        matches_opts=dict(),
        svg_opts=dict(),
        append=False, quiet=False,
        log='log.txt', loglevel='info', redirect_warnings=False):
    '''
    Create Gemoic Comparison Web
    
    Necessary Arguments:
        
        Name               Type     Description
        
        genome_array       list     list of paths to genomes                

    Standard Options:
    
        reference_genome   str      path to reference to order contigs
                                    by (required if reorder=True)
        reference_genome   list     list of reference genomes (matched
                                    against genome_array)
        working_directory  str      path to scratch space for program
                                    to write files
        out_file           str      path to output SVG file
        include_referecne  bool     include referecne genome in
                                    resulting map                            
        reorder            bool     reorder contigs against reference
        palette            str      color palette to use from 
                                    colorcet.palette, e.g. "bgy"
        palette            list     custom palette with list of colors
                                    to use e.g. ["#000000", "red"]. See
                                    svgwrite for colors accepted.
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
                                    
    Basic Geometry Gptions
    
        size               int      size of hive plot in px
        x                  int      x origin, default = 0
        y                  int      y origin, dafault = 0
        rotation           float    value in degrees to rotate all axes
                                    by (clockwise direction)
                                    
    Mulit-panel Geometry Options:
        
        append             bool     append output to existing SVG,
                                    default to append is out_file
        append             str      file to append
        width              int      width of final SVG in px
        height             int      height of fnial SVG in px
        
    Advanced Geometry Options:
    
        dummy_axes         int      number of spacer axes to insert
                                    (useful for pairwise comparisons)
        bezier_max_n       int      max number of genomes before
                                    straight lines are used instead of
                                    bezier curves (set to 0 for always
                                    straight)
        source_angle       float    source angle for bezier curve
                                    (radians)
        target_angle       float    target angle for bezier curve
                                    (radians)
    
    Color Palette Options:
    
        palette_usage      float    decimal percent of palette spectrum
                                    to use
        invert             bool     invert colors in palette
        flip               bool     use palette in reverse order
                    
    Advanced Options:

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
    
    # ======== Initialise and validate options =========

    if log:
        log_level = getattr(logging, loglevel.upper(), None)
        if not isinstance(log_level, int):
            raise ValueError('Invalid log level: %s' % loglevel)
        
        logging.basicConfig(
            filename=join(working_directory, log),
            filemode='w',
            format='%(asctime)s %(message)s',
            datefmt='%d/%m/%Y %I:%M:%S %p',
            level=log_level)
        logging.captureWarnings(redirect_warnings)
    
    # Check all files exist before running anything
    if isinstance(reference_genome, (list, tuple)):
        if reorder:
            assert len(reference_genome) >= len(genome_array), \
            'Unequal number of references and query genomes.'
            # do not raise exception if greater as this will not cause
            # the downstream functions to fail
            if len(reference_genome) > len(genome_array):
                warn('Number of reference genomes greater than input genomes.')
            for f in reference_genome:
                if f:
                    assert exists(f), '%s cannot be found.' % f
    elif reorder:
        assert exists(reference_genome), \
        'Reference genome: %s cannot befound.' % reference_genome
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
    
    # Setup color palette config
    if isinstance(palette, (list, tuple)):
        palette_len = len(palette) - 1
    else:
        try:
            palette = cc.palette[palette]
        except KeyError:
            raise KeyError('Invalid color palette selection: %s.' % palette)
        palette_len = len(palette) - 1
    if invert:
        for i, c in enumerate(palette):
            c_int = int(c[1:], 16)
            n_int = 0xFFFFFF - c_int
            palette[i] = ('#%6x' % n_int).replace(' ', '0').upper()
    if flip:
        palette = palette[::-1]
            
    # ========= Setup and Check Geometry =========
    
    # global geometry
    if width is None:
        width = size
    if height is None:
        height = size
    if viewBox is None:
        viewBox = (0, 0, width, height)
    
    axis_warn = 'Some or all of the axes may not be visible.'
    radii_warn = 'Radii total greater than 100 %. '
    
    # inner-shell radius
    pre_r = int(inner_radius * size / 200)
    # outer-shell radius
    r = int(outer_radius * size / 200) + pre_r
    
    if border_offset:
        # outer-shell boarder offset
        offset = int(border_offset * size / 200)
    else:
        offset = size / 2 - r
        if offset < 0:
            warn(radii_warn + axis_warn, GeometryWarning)
    
    if r + offset > size / 2:
        warn(radii_warn + axis_warn, GeometryWarning)
    
    # origin    
    ox = x + r + offset
    oy = y + r + offset
    ori_warn = 'The origin is outside of the view area. '
    if ox < 0 or ox > width or oy < 0 or oy > height:
        warn(ori_warn, GeometryWarning)
    
    # this may be done on purpose, filter warning
    plot_warn = 'The plot size is greater than the view area. '
    if width < size or height < size:
        warn(plot_warn + axis_warn, GeometryWarning)
    
    
    # ========= Re-index Contigs =========
    
    # set label names
    if not label_names:
        names = [splitext(basename(g))[0] for g in genome_array]
    else:
        assert len(label_names) >= len(genome_array), \
        'Number of labels does not match number of genomes.'
        names = label_names
        if len(label_names) > len(genome_array):
            warn('There are more labels than axes.')
    
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
    # line filtering is also done here
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
    
    # ========= Define Axes =========

    logging.info('Creating axes.')
    positions = []
    if genome_array:
        n_axes = len(genome_array) + dummy_axes
    else:
        n_axes = len(label_names) + dummy_axes
    theta = (2 * np.pi) / n_axes
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
    
    # ========= Define hive plot =========

    # save previous SVG in memory, if in append mode
    if append:
        if isinstance(append, bool):
            template = su.transform.fromfile(out_file)
        else:
            template = su.transform.fromfile(append)
    
    hive = Hiveplot(
        out_file,
        size=('%dpx' % width, '%dpx' % height),
        viewBox='%d %d %d %d' % viewBox,
        **svg_opts)
    
    # link the axes to the plot
    for i, (a, b) in enumerate(positions):
        if i < n_axes - dummy_axes:
            hive.axes.append(Axis(start=a, end=b, **axes_opts))
    
    # ========= Connect up the Axes =========
    
    if 'pid_cov' in matches_opts:
        pid_cov = matches_opts['pid_cov']
    else:
        pid_cov = 90.0
    logging.info('Connecting axes with cutoff of %.1f %%' % pid_cov)
    n_paths = 0
    curved = False if len(genome_array) > bezier_max_n else True

    # default connection angles = 45 deg
    if source_angle is None:
        source_angle = np.pi / 4
    if target_angle is None:
        target_angle = -np.pi / 4

    for i, match in enumerate(matches):
        if i or not dummy_axes:
            for node_index, (pos1, ctg1, pos2, ctg2, pid) in enumerate(match):
                # calculate Node offsets
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
                # convert pid into palette index
                pid = int(round(palette_len * palette_usage * (pid - pid_cov) / (100 - pid_cov)))
                n_paths += 1
                hive.connect(
                    hive.axes[i-1], node_index, source_angle,
                    hive.axes[i], node_index, target_angle,
                    curved=curved,
                    stroke=palette[pid],
                    **connection_opts)
            
            # clear node arrays as they are not needed anymore
            # should save some memory in large graphs
            hive.axes[i-1].nodes = []
            hive.axes[i].nodes = []
    
    # ========= Add Labels to Axes =========

    if not custom_font:
        font = 'font-size:%dpx; font-family:%s' % (font_size, font_family)
    else:
        font = custom_font
    if add_labels:
        for i, ax in enumerate(hive.axes):
            ax.add_node(
                Node(names[i], draw_label=True),
                offset=label_offset, font=font)
    
    # ========= Save and write svg =========

    msg = 'Total number of connections: %d' % n_paths
    logging.info(msg)
    if not quiet:
        print(msg)
        if append:
            print('Saving and merging files...')
        else:
            print('Saving %s...' % out_file)
    logging.info('Output file: %s' % out_file)
    hive.save()
    if append:
        template.append(su.transform.fromfile(out_file))
        template.save(out_file)
    else:
        template = su.transform.fromfile(out_file)
    if not quiet:
        print('Finished.')
    return template