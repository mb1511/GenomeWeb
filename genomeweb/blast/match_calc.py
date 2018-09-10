'''
Created on 15 May 2018

@author: mb1511
@organization: University of Bristol
@contact: mb1511@bristol.ac.uk
@summary:
'''
from __future__ import print_function
from builtins import range

import logging
import re

from genomeweb.blast import local_blast
from os.path import join

def get_matches(
        n1, n2, n1_path='', n2_path='', chunk=2000, step=5000,
        short_ctgs=False, shortck=30, shortsp=100, wd='', l_max_hsps=4,
        quiet=False, **kw):
    
    # check BLAST can be found - will raise error if not
    local_blast.check_blast()

    full = 0
    short = False
    
    with open(join(wd, 'nuc_1.fna'), 'w') as s, \
         open(join(wd, 'nuc_short.fna'), 'w') as d:
        
        logging.info('Creating query files from %s.' % n1_path)
        
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
    
        if full:
            logging.info('%d bp queries: %s' % (chunk, s.name))
        if short:
            logging.info('%d bp queries: %s' % (shortck, d.name))
        
    with open(join(wd, 'nuc_2.fna'), 'w') as s:
        
        logging.info('Reformatting search file %s.' % n2_path)
        
        for i, ctg in enumerate(n2):
            name = re.sub('[\n\r]', '', ctg.name) + '_[ctg=%d]\n' % i
            s.write(name)
            s.write(ctg.seq)
            s.write('\n')
        
        logging.info('Out file:  %s.' % s.name)
    
    # set blast defaults - can be overwritten in **kw
    short_defaults = dict(
        db_path=join(wd, 'db.fna'),
        mr=0, mev=10000,
        join_hsps=False,
        b_type='blastn',
        r_a=True,
        num_threads=4,
        quiet=quiet)
    short_defaults.update(kw)
    
    # override task for short queries
    if short_defaults['b_type'] == 'blastn':
        short_defaults.update(task='blastn-short')
    elif short_defaults['b_type'] == 'blastp':
        short_defaults.update(task='blastp-short')
    
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
        max_hsps=l_max_hsps,
        quiet=quiet)
    
    defaults.update(kw)
    long_defaults.update(defaults)
    
    if short:
    
        local_blast.run(
            join(wd, 'nuc_2.fna'),
            query=join(wd, 'nuc_short.fna'),
            out=join(wd, 'temp_blast_1.xml'),
            outfmt='details',
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
                blast_run=False, make_db=False,
                quiet=quiet)
    
    else:
        itrs = local_blast.run(
            db=join(wd, 'nuc_2.fna'),
            query=join(wd, 'nuc_1.fna'),
            out=join(wd, 'temp_blast.xml'),
            outfmt='details',
            **defaults)
    logging.info('BLAST Search Complete.')
    return itrs