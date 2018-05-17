# Genome Web

## Description

A Python script to create Genomic Webs based on percentage identity

## Installation

Download and install BLAST from [NCBI site](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
(tested with blast-2.7.1+)

If not already installed, run:
	
	pip install numpy colorcet svgwrite

Run:

	./setup.py install

## Usage

	import genomeweb
	
	if __name__ == '__main__':
		genomeweb.create_web(genome_list, reference_genome, **options)

If BLAST cannot be found on the system path, set it explicitly before running:

	genomeweb.blast.BLAST_PATH = path/to/NCBI/blast-x.x.x+/bin
		
Necessary Arguments:

	Name               Type     Description                             Default
	
	genome_array       list     list of paths to genomes				
	reference_genome   str      path to reference to order contigs by
	                            (not required if reorder=False)

Standard Options:

	working_directory  str      path to scratch space for program to
	                            write files
	out_file           str      path to output SVG file                 'default.svg'
	include_referecne  bool     include referecne genome in resulting   False
	                            map							
	reorder            bool     reorder contigs against reference       True
	palette            str      color palette to use from cc.palette    'bgy'
	add_labels         bool     add labels to axes                      True
	label_offset       float    label offset                            1.1
	inner_radius       int      inner-shell radius                      100
	outer_radius       int      distance from inner shell to outer-     300
	                            shell bounds
	border_offset      int      distance from border (increase to       80
	                            correctly display longer genome names)
	matches_opts       dict     kwargs to parse to find matches 
	                            function, such as %id cutoff value
	                            (pid_cov) and minimum hit length
	                            (hit_len) for use in filtering hits
	                            
Advanced Options:

	palette_usage      float    decimal percent of palette spectrum     1.0
	                            to use
	bezier_max_n       int      max number of genomes before straight   4
	                            lines are used instead of bezier
	                            curves (set to 0 for always striaght)
	x_scaling          str      SVG x scaling factor                    '200%'
	y_scaling          str      SVG y scaling factor                    '200%'
	connection_opts    dict     dictionary of line options
	axes_opts          dict     dictionary of axes options
	reorder_opts       dict     kwargs to parse to reorder function
	
	


## Other Notes

The pyveplot2 module is a heavily modified version of [pyveplot](https://github.com/rgarcia-herrera/pyveplot)


