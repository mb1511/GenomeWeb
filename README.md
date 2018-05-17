# Genome Web

# Description

A Python script to create Genomic Webs based on percentage identity

# Installation

Download and build the BLAST package from NCBI ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

If not already installed, run:
	
	pip install numpy, colorcet

Run:

	./setup.py install

# Usage

	import genomeweb
	
	if __name__ == '__main__':
		genomeweb.create_web(genome_list, reference_genome, **options)

If BLAST cannot be found on the system path, set it explicitly before running:

	genomeweb.blast.BLSAT_PATH = path/to/NCBI/blast+x.x.x/bin
		
Necessary Args:

		
Other args:


	
See examples for more information.

# Other Notes:


