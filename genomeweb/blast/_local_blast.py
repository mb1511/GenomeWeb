'''
Created on 11 Jun 2016

@author: Matt

@todo: fix db_path for no db selected -change to default
'''
from __future__ import print_function 

import subprocess as sp
from os.path import join
import logging, warnings

from genomeweb.blast import read

BLAST_PATH = ''

DB_TPYE = {
    'blastp': 'prot',
    'blastx': 'prot',
    'blastn': 'nucl',
    'tblastn': 'nucl',
    'tblastx': 'nucl'
}

def check_blast():
    global BLAST_PATH
    # try a call to blastp
    logging.debug('Checking BLSAT install.')
    try:
        cmd = [
            'blastp' if not BLAST_PATH else join(BLAST_PATH, 'blastp'),
            '-version']
        logging.debug(' '.join(cmd))
        o, e = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE).communicate()
        logging.debug(o)
        logging.debug(e)
    except OSError:
        logging.debug('BLAST not found.')
        raise OSError('BLAST executables not found. Please set genomeweb.blast.BLAST_PATH to path/to/NCBI/blast-x.x.x+/bin')
    logging.debug('BLAST found.')
    return True

def run(
        db=None, db_path='', query='', subject='', out='-',
        save_subject_db=False, make_db=True, db_type=None, b_type='blastp',
        outfmt='clustal', mev=0.001, mr=0, join_hsps=False, 
        blast_run=True, r_a=False, quiet=False, *args, **kwargs):
    '''
    Local BLAST
    
    db      -    raw .faa/fna proteome/genome file path - formats blast
                 db from file if `make_db=True`
    db_path -    path to blast database file (minus extension) - will be
                 created if `make_db=True`
    
    query   -    path to query sequence; must be fasta file, can have 
                 multiple entries
    subject -    genome/proteome file to search, does not save a database
    
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
    global BLAST_PATH, DB_TPYE
        
    if BLAST_PATH:
        makeblastdb = join(BLAST_PATH, 'makeblastdb')
        blast = join(BLAST_PATH, b_type)
    else:
        makeblastdb = 'makeblastdb'
        blast = b_type
    
    if blast_run:
        assert db or subject or (db_path and not make_db), 'No subject sequence or database file given.'
        assert query, 'No query sequence file provided.'

    if db_type is None:
        # try to auto assign database molecule type
        try:
            db_type = DB_TPYE[b_type]
        except KeyError:
            # in case of other non normal BLAST options
            msg = 'Automatic database molecule assignment not configured for: "{}".\
            Falling back to db_type="nucl".'.format(b_type)
            warnings.warn(msg)
            logging.info(msg)
            db_type = 'nucl'
    
    if make_db or save_subject_db:
        if not db and subject and save_subject_db:
            db = subject
        if db:
            # copy file to new location
            with open(db, 'r') as s:
                with open(db_path, 'w') as q:
                    q.write(s.read())
                    logging.info('Copying %s to %s.' % (db, db_path))
            
            # build database using makeblastdb
            cmd = [makeblastdb, '-in', db_path, '-parse_seqids','-dbtype', db_type]
            mdp = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
            o, e = mdp.communicate()
            if not quiet:
                print(o)
                print(e)
            logging.debug('--- BEGIN MAKEBLASTDB ---')
            logging.info(' '.join(cmd))
            logging.debug('StdOut:\n%s\nStdErr\n%s' % (o, e))
            logging.debug('--- END MAKEBLASTDB ---')
    
    # append any additional arguments to command
    add_args = []
    for arg in args:
        add_args.append('-' + str(arg))
    for key in kwargs.keys():
        if isinstance(kwargs[key], bool):
            if kwargs[key]:
                add_args.append('-' + str(key))
        else:
            add_args.extend(['-' + str(key), str(kwargs[key])])
    
    # perform BLAST search
    if blast_run:
        # will just return what ever is in the temp_blast.xml if False  
        if not quiet:      
            print('Running BLAST: %s' % blast)

        cmd = [blast, '-query', query, '-out', out, '-outfmt', '5'] + add_args
        if db_path:
            cmd.extend(['-db', db_path])
        elif subject:
            cmd.extend(['-subject', subject])
        else:
            raise RuntimeError('No database or subject sequence file path is given.')

        _b = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT)
        o, e = _b.communicate()
        logging.debug('--- BEGIN BLAST ---')
        logging.info(' '.join(cmd))
        logging.debug('StdOut:\n%s\nStdErr\n%s' % (o, e))
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