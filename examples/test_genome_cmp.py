'''
Script for generating examples
'''
import glob
import os
from os.path import basename
import genomeweb as gw

# check for blast bin
if 'blast-' not in os.environ['PATH']:
    gw.blast.BLAST_PATH = 'path/to/blast/bin'

# For non-disclosure reasons, the genomes used here have not been added
# on GitHub yet, but will be in the near future.

# A list of 20 genomes with names R______.fna
files = [
    f for f in glob.glob('genomes/*.fna') if basename(f).startswith('R')]

# reference can be anything - but this a type strain in the same genus
# as the others
reference = 'genomes/reference.fna'

# tri-spoke, no re-indexing
gw.create_web(
    files[:3],
    reference_genome=reference,
    working_directory='./scratch',
    out_file='unordered_tri.svg',
    reorder=False,
    palette_usage=0.8)

# tri-spoke, with re-indexing
gw.create_web(
    files[:3],
    reference_genome=reference,
    working_directory='./scratch',
    out_file='tri.svg',
    palette_usage=0.8)

# 20-spoke, with re-indexing
gw.create_web(
    files,
    reference_genome=reference,
    working_directory='./scratch',
    out_file='unordered_tri.svg')

