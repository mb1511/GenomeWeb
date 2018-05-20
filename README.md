# Genome Web

## Description

A Python script to create Genomic Webs based on percentage identity.

## Installation

Requires Python 2.7. Working on a 3.6-compatible version.

Download and install BLAST from <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+> (tested with blast-2.7.1+)

If using an installer, it should automatically prepend the system path.

If the following libraries are not already installed, run the following command from a Terminal/Command Prompt:
	
	pip install numpy colorcet svgwrite svgutils

(pip should be installed by default with Python)

Install the scripts:

* Download and unzip this repository.
* Navigate to the download direcotry and open a Terminal/Command Prompt/PowerShell here (Shift + Right-click on Windows gives the option for this)
* Run the install script:

On UNIX-based systems:

	./setup.py install
	
On any platform:

	python setup.py install

## Usage

	import genomeweb
	
	genomeweb.create_web([
			'genome1.fna',
			'genome2.fna',
			'genome3.fna'],		# can be a list of arbitrary length
		'reference_genome.fna',	# can be any genome (only necessary when re-indexing)
		**options)			# see list of available options below

If BLAST cannot be found on the system path, set it explicitly before running:

	genomeweb.blast.BLAST_PATH = path/to/NCBI/blast-x.x.x+/bin

An `svgutils.transform.SVGFigure` is returned for any annotations to be made post plot completion. See [svgutils](https://github.com/btel/svg_utils) for full documentation.
		
### Necessary Arguments
        
        Name               Type     Description
        
        genome_array       list     list of paths to genomes                

### Standard Options
    
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
                                    
### Basic Geometry Options
    
        size               int      size of hive plot in px
        x                  int      x origin, default = 0
        y                  int      y origin, dafault = 0
        rotation           float    value in degrees to rotate all axes
                                    by (clockwise direction)
                                    
### Mulit-panel Geometry Options
        
        append             bool     append output to existing SVG,
                                    default to append is out_file
        append             str      file to append
        width              int      width of final SVG in px
        height             int      height of fnial SVG in px
                        
### Advanced Options

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
    
### Additional Options (likely do not need to be altered)
    
        border_offset      int      distance from border, default
                                    calculated automatically
        viewBox            tuple    viewBox for SVG, default
                                    calculated automatically  
	                            
### Example for editing match searching and filtering options

	create_web(
		files, ref,
		matches_opts=dict(
			chunk=3000,		# query length of 3000 bp
			step=2000,		# step length of 2000 bp between queries
			pid_cov=85,		# filter out hits below 85 % identity
			hit_len=1200))		# filter out hits less than 1200 bp

CAUTION: keep filtering resonably high or drastically increase step length to increase speed, clarity and reduce file sizes


If adding connection or axes options, the defaults will all be overwritten. 

### Default connection options

	connection_opts=dict(stroke_width='0.34', stroke_opacity='0.4')
	
### Default axes options

	axes_opts=dict(stroke='black', stroke_width='1')
	
### Example for reorder_opts for increase reordering speed (similar set of options to match finding)

	create_web(
		files, ref,
		reorder_opts=dict(
			chunk=2000,		# query length of 2000 bp
			step=10000,		# step length of 10000 bp between queries
			short_ctgs=False))	# throws away contigs less than chunk size

### Example for reorder_opts but retaining short contigs

	create_web(
		files, ref,
		reorder_opts=dict(
			chunk=2000,		# query length of 2000 bp
			step=10000,		# step length of 10000 bp between queries
			shortck=30,		# short contig query of 30 bp
			shortsp=200))		# short contig step of 200 bp

NOTE: if a hit for any position within a contig cannot be found, the whole contig will be discarded. Contigs will also be discarded if their length is below the minimum query size and when short contig searching is turned off (on by default).

If the query genome is very distant to the reference (i.e. no regions of homology - very unlikely), contig disposal can be completely turned off if reorder is set to False to prevent any sorting, but the resulting plot will have very few or no connections at all.

### Adding Annotations to the Plot
	
	import genomeweb as gw
	import svgutils.transform as sg
	fig = gw.create_web(genome_list, out_file='figure.svg', **options)
	# e.g. add an "A" to the top left of the figure.
	fig.append(sg.TextElement(0, 20, 'A', size=12))
	fig.save('figure.svg')

### Creating Multi-panel Figures

	# "width" is the final size with both panels (1000 for two adjacent 500px panels)
	# "x" is the x anchor point for the panel
	# "height" and "y" can be used in the same way
	options = dict(out_file='figure.svg', size=500, width=1000, reorder=False)
	# add first panel at 0,0 
	gw.create_web(genome_list1, **options)
	# add second panel at 500,0
	gw.create_web(genome_list2, x=500, append=True, **options)
	
	

## Citations

If used for academic publications, please cite the following publications:

*A Publication Comming Soon Here*. In the meantime, please cite the GitHub page.

Camacho C., *et al* (2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/20003500?dopt=Citation)

## Other Notes

The pyveplot2 module is a heavily modified version of [pyveplot](https://github.com/rgarcia-herrera/pyveplot)

## Examples

See example [script](https://github.com/mb1511/GenomeWeb/blob/master/examples/test_genome_cmp.py) for generation of the following plots:

### Plot with no contig re-indexing:

![unordered](https://github.com/mb1511/GenomeWeb/blob/master/examples/unordered_tri.svg)

### Plot post contig re-indexing:

![three spoke](https://github.com/mb1511/GenomeWeb/blob/master/examples/tri.svg)

### 20 Just for fun:

![twenty](https://github.com/mb1511/GenomeWeb/blob/master/examples/twenty.svg)

### Multi-panel Example.
See [multi_panel.py](https://github.com/mb1511/GenomeWeb/blob/master/examples/multi_panel.py) for source.

![multi](https://github.com/mb1511/GenomeWeb/blob/master/examples/multi_panel.svg)
