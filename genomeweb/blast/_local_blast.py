'''
Created on 11 Jun 2016

@author: Matt

@todo: fix db_path for no db selected -change to default
'''
from __future__ import print_function 

import subprocess as sp
import os
from genomeweb.blast import read

BLAST_PATH = ''

if 'blast-' not in os.environ['PATH']:
    raise OSError('BLAST executables not found. Please set genomeweb.blast.BLAST_PATH to path/to/NCBI/blast-x.x.x+/bin')
    #BLAST_PATH = 'path/to/NCBI/blast-2.7.1+/bin/'

def run(db=None, db_path='', 
        query = '',
        out = '-', 
        mev=0.001, mr=0, join_hsps=False, b_type='blastp', 
        make_db=True, db_type=None, outfmt='clustal',
        blast_run=True, r_a=False, quiet=False, *args, **kwargs):
    u'''
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
        
        # build database using makeblastdb.exe        
        mdp = sp.Popen([makeblastdb, '-in', db_path,
                       '-parse_seqids','-dbtype', db_type],
                       stdout=sp.PIPE, stderr=sp.PIPE)
        o, e = mdp.communicate()
        if not quiet:
            print(o)
            print(e)
    
    add_args = []
    for arg in args:
        add_args.append('-' + str(arg))
    for key in kwargs.keys():
        add_args.extend(['-' + str(key), str(kwargs[key])])
    
    # perform BLAST using blast exe
    if blast_run:
        # will just return what ever is in the temp_blast.xml if False  
        if not quiet:      
            print('Running BLAST: %s' % blast)
        _o = sp.Popen([blast, '-query', query,
                       '-db', db_path, '-out', out, '-outfmt', '5']
                      + add_args, stdout=sp.PIPE, stderr=sp.STDOUT)
        o = _o.communicate()[0]
        #print(o)
        #print(e)
    
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
        return 'No match found.'
    
    return output






