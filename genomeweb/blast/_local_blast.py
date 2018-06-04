'''
Created on 11 Jun 2016

@author: Matt

@todo: fix db_path for no db selected -change to default
'''
from __future__ import print_function 

import subprocess as sp
import os
import logging

from genomeweb.blast import read

BLAST_PATH = ''

# can override blast check if required
OK = False


def run(db=None, db_path='', 
        query = '',
        out = '-', 
        mev=0.001, mr=0, join_hsps=False, b_type='blastp', 
        make_db=True, db_type=None, outfmt='clustal',
        blast_run=True, r_a=False, quiet=False, *args, **kwargs):
    '''
    Local BLAST
    
    db      -    raw *.faa/fna proteome file path; formats blast db from 
                 file if make_db is True
    db_path -    path to blast data file
    
    query   -    path to query sequence; must be fasta file, can have 
                 multiple entries
    
    mev     -    max E-val; set to 1 to retrun all
    mr      -    max alignments to return; default = 0 -> returns all
    
    make_db -    set to false to use current stored temp blast db
    db_type -    database type to create:
                    - None (default); works out automatically
                    - prot
                    - nuc
                    
    outfmt  -    output format:
                    - clustal (default); standard clustal alignment 
                      format
                    - fasta; FASTA format
                    - details; returns details only in tuple with no 
                      formatting
    
    join_hsps -  True/False; join contiguous hsps in BLAST hits;
                 if False -> only returns highest scoring hsp
    
    b_type  -    BLAST algorithm to use:
                    - blastp (default)
                    - blastn
                    - blastx
                    - tblastn
                    - tblastx
                    - psiblast
                    - rpsblast
                    - rpsblastn
                    - deltablast
                    
    blast_run -  Set to False if already using premade xml output (skips BLAST execution)
    
    r_a      -   Return all hsps as lists
                    
    *args    -   additional arguments to parse to blast application 
                 e.g: h -> adds argument [..., '-h']
    **kwargs -   additional arguments to parse to blast application 
                 e.g: arg = value -> adds argument [..., '-arg', 'value'] 
    '''
    global OK, BLAST_PATH
    
    if not OK:
        # check path first
        if 'blast-' not in os.environ['PATH']:
            # try a call to blastp
            try:
                sp.call(['blastp', '-version'], stdout=sp.PIPE, stderr=sp.PIPE)
            except OSError:
                raise OSError('BLAST executables not found. Please set genomeweb.blast.BLAST_PATH to path/to/NCBI/blast-x.x.x+/bin')
        OK = True
    
    makeblastdb = os.path.join(BLAST_PATH, 'makeblastdb')
    blast = BLAST_PATH + b_type
    
    if db_type is None:
        if b_type.startswith('t') or 'n' in b_type:
            db_type = 'nucl'
        else:
            db_type = 'prot'
    
    if make_db:
        # copy query.faa file to temp.faa
        with open(db, 'r') as s:
            with open(db_path, 'w') as q:
                q.write(s.read())
                logging.info('Copying %s to %s.' % (db, db_path))
        
        # build database using makeblastdb.exe
        cmd = [makeblastdb, '-in', db_path, '-parse_seqids','-dbtype', db_type]
        mdp = sp.Popen(cmd,
                       stdout=sp.PIPE, stderr=sp.PIPE)
        o, e = mdp.communicate()
        if not quiet:
            print(str(o))
            print(str(e))
        logging.debug('--- BEGIN MAKEBLASTDB ---')
        logging.info(' '.join(cmd))
        logging.debug('StdOut:\n%s\nStdErr\n%s' % (str(o), str(e)))
        logging.debug('--- END MAKEBLASTDB ---')
    
    add_args = []
    for arg in args:
        add_args.append('-' + str(arg))
    for key in kwargs.keys():
        if isinstance(kwargs[key], bool):
            if kwargs[key]:
                add_args.append('-' + str(key))
        else:
            add_args.extend(['-' + str(key), str(kwargs[key])])
    
    # perform BLAST using blast exe
    if blast_run:
        # will just return what ever is in the temp_blast.xml if False  
        if not quiet:      
            print('Running BLAST: %s' % blast)
        cmd = [
            blast, '-query', query, '-db', db_path, '-out', out, '-outfmt', '5'] + add_args
        _b = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT)
        o, e = _b.communicate()
        logging.debug('--- BEGIN BLAST ---')
        logging.info(' '.join(cmd))
        logging.debug('StdOut:\n%s\nStdErr\n%s' % (str(o), str(e)))
        logging.debug('--- END BLAST ---')
        
    
    if out != '-':
        rd = False
    else:
        rd = True
        out = o
    
    try:
        if outfmt == 'clustal':
            output = read.clustal(out, raw_data=rd, num_ret=mr, expect=mev)
        elif outfmt == 'fasta':
            output = read.fasta(out, raw_data=rd, num_ret=mr, expect=mev)
        elif outfmt == 'details':
            # returns list [ query0[ hit0, ... ], ... ]
            output = read.details(
                out, raw_data=rd, num_ret=mr,
                join_hsps=join_hsps, expect=mev, r_a=r_a)
    except StopIteration:
        # null record
        logging.info('No match found in BLAST search.')
        return 'No match found.'
    
    return output






