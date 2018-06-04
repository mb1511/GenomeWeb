# Genome Web

## Description

A Python script to create Genomic Webs based on percentage identity.

## Contents

* [Installation](https://github.com/mb1511/GenomeWeb#installation)
* [Usage](https://github.com/mb1511/GenomeWeb#usage)
	* [Necessary Arguments](https://github.com/mb1511/GenomeWeb#necessary-arguments)
	* [Standard Options](https://github.com/mb1511/GenomeWeb#standard-options)
	* [Basic Geometry Options](https://github.com/mb1511/GenomeWeb#basic-geometry-options)
	* [Multi-panel Geometry Options](https://github.com/mb1511/GenomeWeb#multi-panel-geometry-options)
	* [Advanced Geometry Options](https://github.com/mb1511/GenomeWeb#advanced-geometry-options)
	* [Color Palette Options](https://github.com/mb1511/GenomeWeb#color-palette-options)
	* [Advanced Options](https://github.com/mb1511/GenomeWeb#advanced-options)
	* [BLAST+ Search Options](https://github.com/mb1511/GenomeWeb#blast-search-options)
	* [Additional Options](https://github.com/mb1511/GenomeWeb#additional-options-likely-do-not-need-to-be-altered)
	* [Adding Annotations](https://github.com/mb1511/GenomeWeb#adding-annotations-to-the-plot)
	* [Creating Multi-panel Figures](https://github.com/mb1511/GenomeWeb#creating-multi-panel-figures)

* [Citations](https://github.com/mb1511/GenomeWeb#citations)
* [Other Notes](https://github.com/mb1511/GenomeWeb#other-notes)
* [Examples](https://github.com/mb1511/GenomeWeb#other-notes)

## Installation

Requires Python 2.7 or 3.x. Originally built in 2.7, but have added 3.x compatibility. Tested with 2.7 and 3.6.

Download and install BLAST from <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+> (tested with blast-2.7.1+)

If using an installer, it should automatically prepend the system path.

If the following libraries are not already installed, run the following command from a Terminal/Command Prompt:
	
	pip install numpy colorcet svgwrite svgutils
	
I using Python 2.7.x, the following may also need to be installed:

	pip install future

(pip should be installed by default with Python)

Install the scripts:

* Download and unzip this repository.
* Navigate to the download directory and open a Terminal/Command Prompt/PowerShell here (Shift + Right-click on Windows gives the option for this)
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
		'reference_genome.fna',		# can be any genome (only necessary when re-indexing)
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
                                    
### Basic Geometry Options
    
        size               int      size of hive plot in px
        x                  int      x origin, default = 0
        y                  int      y origin, dafault = 0
        rotation           float    value in degrees to rotate all axes
                                    by (clockwise direction)
                                    
### Multi-panel Geometry Options
        
        append             bool     append output to existing SVG,
                                    default to append is out_file
        append             str      file to append
        width              int      width of final SVG in px
        height             int      height of fnial SVG in px

### Advanced Geometry Options
    
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

### Color Palette Options
    
        palette_usage      float    decimal percent of palette spectrum
                                    to use
        invert             bool     invert colors in palette
        flip               bool     use palette in reverse order

### Advanced Options

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
	
### BLAST+ Search Options

For advanced users who wish to alter the default BLAST search options, command line options can be parsed directly to BLAST via key-word arguments in the `matches_opts` and `reorder_opts` dictionaries. For example:
	
	gw.create_web(
		files, ref,
		matches_opts=dict(
			word_size=15))		# appends "-word_size 15" to blast command

To run megablast, dicontiguous megablast or blastn, use option `task=` either `megablast`, `dc-megablast` or `blastn`. Default is `blastn`.

Below is a list of arguments that can be parsed to blastn. Do not alter `query`, `db`, `out`, `outfmt` or `subject`. The database path (`db`) can however be modified using the `db_path="path"` argument with `make_db=False`. This will then not automatically create a blast database, but will use the one specified. This is probably useless for the `matches_opts` as currenty only one custom database can be specified. It is intended for use in the `reorder_opts`, though this will also limit the user to one reference genome.

The number of threads to run blast on can be set using the `num_threads=int` option. Default is 4, though it will likely not fully leverage these over the short runtime.

To edit `-max_hsps` for full-length searches (not short-contig searching) use the `l_max_hsps=` option; default is 4. To edit this for short-contigs, just add `max_hsps=val`; default is not to include this and return all possible hsps.

For positional arguments with no value, e.g. `-lcase_masking` use `lcase_masking=True` to add the argument. In this instance, lower-case basepairs (typically low confidence sequences) are masked in the query and subject sequences .

Options for blastn:

	  blastn [-h] [-help] [-import_search_strategy filename]
	    [-export_search_strategy filename] [-task task_name] [-db database_name]
	    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
	    [-negative_gilist filename] [-entrez_query entrez_query]
	    [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
	    [-subject subject_input_file] [-subject_loc range] [-query input_file]
	    [-out output_file] [-evalue evalue] [-word_size int_value]
	    [-gapopen open_penalty] [-gapextend extend_penalty]
	    [-perc_identity float_value] [-qcov_hsp_perc float_value]
	    [-max_hsps int_value] [-xdrop_ungap float_value] [-xdrop_gap float_value]
	    [-xdrop_gap_final float_value] [-searchsp int_value]
	    [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy]
	    [-min_raw_gapped_score int_value] [-template_type type]
	    [-template_length int_value] [-dust DUST_options]
	    [-filtering_db filtering_database]
	    [-window_masker_taxid window_masker_taxid]
	    [-window_masker_db window_masker_db] [-soft_masking soft_masking]
	    [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
	    [-best_hit_score_edge float_value] [-window_size int_value]
	    [-off_diagonal_range int_value] [-use_index boolean] [-index_name string]
	    [-lcase_masking] [-query_loc range] [-strand strand] [-parse_deflines]
	    [-outfmt format] [-show_gis] [-num_descriptions int_value]
	    [-num_alignments int_value] [-line_length line_length] [-html]
	    [-max_target_seqs num_sequences] [-num_threads int_value] [-remote]
	    [-version]

See [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK279668/) for BLAST+ manual and [here](https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a) for all command line options explained. Use only if there is a good reason, otherwise defaults should work fine.

`-perc_identity val` (in BLAST) filtering will effectively be overriden by the `pid_cov` value (in match filtering post-blast) if `pid_cov > perc_identity`.

Theoretically, the blast protocol can be altered to use a different blast type by adding the `b_type=` argument: `"blastp"`, `"psiblast"`, `"deltablast"` or `"tblastx"` for protein sequences also adding the `db_type="prot"` argument to correctly make the databases. Cross-molecule searches are not supported: `"tblastn"` or `"blastx"`. CAUTION: this is not tested.

## Citations

If used for academic publications, please cite the following publications:

*A Publication Comming Soon Here*. In the meantime, please cite the GitHub page.

Camacho C., *et al* (2008) BLAST+: architecture and applications. BMC Bioinformatics 10:421. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/20003500?dopt=Citation)

## Other Notes

The pyveplot2 module is a heavily modified version of [pyveplot](https://github.com/rgarcia-herrera/pyveplot)

## Examples

See example [script](https://github.com/mb1511/GenomeWeb/blob/master/examples/test_genome_cmp.py) for generation of the following plots:

### Plot with no contig re-indexing

![unordered](https://github.com/mb1511/GenomeWeb/blob/master/examples/unordered_tri.svg)

### Plot post contig re-indexing

![three spoke](https://github.com/mb1511/GenomeWeb/blob/master/examples/tri.svg)

### 20 Just for fun

![twenty](https://github.com/mb1511/GenomeWeb/blob/master/examples/twenty.svg)

### Multi-panel Example

See [multi_panel.py](https://github.com/mb1511/GenomeWeb/blob/master/examples/multi_panel.py) for source.

![multi](https://github.com/mb1511/GenomeWeb/blob/master/examples/multi_panel.svg)

### Pairwise Comparison

See [script](https://github.com/mb1511/GenomeWeb/blob/master/examples/test_double.py) for options to create this.

![pair](https://github.com/mb1511/GenomeWeb/blob/master/examples/double.svg)