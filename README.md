<!-- ABOUT THE PROJECT -->
<a name="readme-top"></a>
## About UPC

## Introduction

This is a Snakemake workflow for evaluating pre-designed diagnostic primers (and probes) against a comprehensive database of genomes (NCBI).

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<a name="installation-top"></a>
## Installation

1. Make sure you have [Mamba](https://github.com/conda-forge/miniforge) and [Snakemake] (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) installed.

2. Clone the repo into your working directory
   ```sh
   git clone https://github.com/sjlow23/upc.git
   cd upc
   ```

<p align="right">(<a href="#installation-top">back to top</a>)</p>


<a name="config-top"></a>
## Configuration file

Modify the configuration file to include the parameters specific for your run. Specifically, the following input files are required:

    1. List of primers (and probes) to be evaluated in tab-delimited format:

        column 1: name of the primer-probe set
        column 2: forward primer sequence (can include IUPAC codes)
        column 3: reverse primer sequence (can include IUPAC codes)
        column 4: probe sequence (optional; can include IUPAC codes)

Place the primer input file into the `input` directory.
The remaining parameters are to be specified in the `config.yaml` file. Details for each parameter as follows:

        `SUBSAMPLE_TARGET` and `SUBSAMPLE_OFFTARGET`: should subsampling be performed for target and off-target database genomes
        `PROBES`: are probes to be evaluated?
        `MIN_GENOME_SIZE`: minimum target genome size for genome to be evaluated
        `DB`: if downloading genomes, should the GenBank or RefSeq database be used?
        `TARGET_SP_TAXID` and `OFFTARGET_SP_TAXID`: if downloading genomes, species taxids to use
        `DOMAIN`: pipeline currently only working for viral genomes
        `ASSEMBLY_LEVEL`: one of "complete", "chromosome", "scaffold", "contig"
        `USE_ASSEMBLY`: should genomes be downloaded from the assembly resource (GCA|GCF)
        
        isPcr parameters:
        `MAX_AMPLICON_SIZE`: maximum product size a primer pair should amplify
        `MIN_PERFECT`: minimum size of perfect match at 3' end of primer (default 15)
        `TILE_SIZE`: the size of match that triggers an alignment (default 11)
        `STEP_SIZE`: spacing between tiles (default 5)
        `MIN_GOOD`: minimum size where there must be 2 matches for each mismatch (default 15)

        Amplicon parameters:
        `MAX_MISMATCHES`: Currently NOT USED

        Probe parameters:
        `MAX_MISMATCHES`: maximum allowed mismatches in probe binding site to be considered "present"



<a name="run-top"></a>
## Run pipeline

   ```sh
   snakemake --cores 16 --use-conda --configfile/config.yaml 
   ```


<a name="output-top"></a>
## Output files

The tables and figures resulting from a UPC run is collated in a HTML document. The details for the tables and figures are listed below. Note that some tables or figures may not be generated if not relevant to dataset.

| Output      	| Details                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         	|
|-------------	|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| Table 1     	| Original list of primers (and probes) provided to be evaluated. Degenerate primer and probe positions are highlighted.                                                                                                                                                                                                                                                                                                                                                                                                                          	|
| Table 2     	| All possible combinations of forward and reverse primers derived from Table 1.                                                                                                                                                                                                                                                                                                                                                                                                                                                                  	|
| Table 3     	| All observed combinations of mismatches observed in the forward primer binding sites, ordered by decreasing prevalence. <br>Mismatches in primer binding sites are highlighted in red. <br>Number of genomes in database refer to original number of target genomes evaluated. <br>Number of genomes amplified refer to genomes successfully amplified _in silico_. <br>Number of genomes with >=1 mismatches refer to count of genomes with at least one mismatch. <br>List of genomes with observed mismatches are listed in the last column. 	|
| Table 4     	| All observed combinations of mismatches observed in the reverse primer binding sites, ordered by decreasing prevalence. <br>Mismatches in primer binding sites are highlighted in blue. Other values are as in Table 3.                                                                                                                                                                                                                                                                                                                         	|
| Table 5     	| All probes derived from Table 1 that were evaluated.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            	|
| Table 6     	| All observed combinations of mismatches observed in the probe binding sites, ordered by decreasing prevalence. <br>Number of genomes with probe binding site refer to number of amplicons successfully amplified that contain the probe binding site. Other values are as in Tables 3 & 4.                                                                                                                                                                                                                                                      	|
| Table 7     	| List of genomes that were not amplified _in silico_ by the primers.                                                                                                                                                                                                                                                                                                                                                                                                                                                                             	|
| Table 8     	| List of genomes with amplicons that did not contain the probe binding site.                                                                                                                                                                                                                                                                                                                                                                                                                                                                     	|
| Figure 1    	| Profiles of primer _in silico_ amplification results.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           	|
| Figure 2    	| Profiles of probe _in silico_ amplification results.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            	|
| Figure 3    	| Prevalence of mismatch combinations observed in forward and reverse primer binding sites. Similar to data in Tables 3 and 4.                                                                                                                                                                                                                                                                                                                                                                                                                    	|
| Figure 4    	| Proportion of genomes with mismatches in the forward and reverse primer binding sites at individual positions.                                                                                                                                                                                                                                                                                                                                                                                                                                  	|
| Figure 5    	| Prevalence of mismatch combinations observed in probe binding sites. Similar to data in Table 6.                                                                                                                                                                                                                                                                                                                                                                                                                                                	|
| Figure 6    	| Proportion of genomes with mismatches in the probe binding sites at individual positions.                                                                                                                                                                                                                                                                                                                                                                                                                                                       	|
| Figure 7    	| Alluvial plot of genome profiles (_in silico_ matches) for primers and probes.                                                                                                                                                                                                                                                                                                                                                                                                                                                                  	|
| Figure 8    	| For genomes where _in silico_ amplification failed, the mismatches observed at the supposed primer binding sites based on multiple sequence alignment with a 'successfully amplified' genome.                                                                                                                                                                                                                                                                                                                                                   	|
| Figure 9    	| For genomes where probe binding sites were not present, the mismatches observed at the supposed probe binding sites based on multiple sequence alignment with a genome containing the probe binding site. <br>Note: genomes with perfect matches to the probes may sometimes be present and still shown in figure. <br>This would be due to a genome that failed _in silico_ amplification but still contained the probe binding site.                                                                                                          	|