configfile: "config/config_denv.yaml"


OUTDIR = config["OUTDIR"]
CPU = config["CPU"]
SUBSAMPLE_TARGET = config["SUBSAMPLE_TARGET"]
SUBSAMPLE_OFFTARGET = config["SUBSAMPLE_OFFTARGET"]
PROBES = config["PROBES"]

PRIMERS = config["PRIMERS"]
MAX_AMPLICON_SIZE = config["MAX_AMPLICON_SIZE"]
MIN_PERFECT = config["MIN_PERFECT"]
TILE_SIZE = config["TILE_SIZE"]
STEP_SIZE = config["STEP_SIZE"]

DOWNLOAD_TARGET = config["DOWNLOAD_TARGET"]
DOWNLOAD_OFFTARGET = config["DOWNLOAD_OFFTARGET"]
DEREP = config["DEREP"]
DATABASE = "database/"
DOMAIN = config["DOMAIN"]
ASSEMBLY_LEVEL = config["ASSEMBLY_LEVEL"]
DB = config["DB"]
USE_ASSEMBLY = config["USE_ASSEMBLY"]
IS_FLU = config["IS_FLU"]

if "USER_TARGET" in config:
	USER_TARGET = config["USER_TARGET"]
if "USER_OFFTARGET" in config:
	USER_OFFTARGET = config["USER_OFFTARGET"]
if "TARGET_SP_TAXID" in config:
	TARGET_SP_TAXID = config["TARGET_SP_TAXID"]
if "OFFTARGET_SP_TAXID" in config:
	OFFTARGET_SP_TAXID = config["OFFTARGET_SP_TAXID"]
if "TARGET_TAXID" in config:
	TARGET_TAXID = config["TARGET_TAXID"]
if "OFFTARGET_TAXID" in config:
	OFFTARGET_TAXID = config["OFFTARGET_TAXID"]

## Need check statement- only one of TARGET_SP_TAXID and TARGET_TAXID can be present
## Need check statement- only one of OFFTARGET_SP_TAXID and OFFTARGET_TAXID can be present

GENOMES_TARGET = OUTDIR + "genomes_target/"
GENOMES_OFFTARGET = OUTDIR + "genomes_offtarget/" 

# Use as target in rule all to get list of target genomes
def get_target_genomes(wildcards):
	target_dir = checkpoints.download_target.get(**wildcards).output[0]
	#global GENOMES_T
	filelist = expand(os.path.join(target_dir, "{genome}.fna"), 
					genome = glob_wildcards(os.path.join(target_dir, "{genome}.fna")).genome)
	# GENOMES_T, = glob_wildcards(os.path.join(target_dir, "{genome}.fna"))
	# filelist = expand(os.path.join(target_dir, "{genome}.fna"), genome = GENOMES_T.genome)
	return filelist

# Use as target in rule all to get list of off-target genomes
def get_offtarget_genomes(wildcards):
	offtarget_dir = checkpoints.download_offtarget.get(**wildcards).output[0]
	filelist = expand(os.path.join(offtarget_dir, "{genome}.fna"), 
					genome = glob_wildcards(os.path.join(offtarget_dir, "{genome}.fna")).genome)
	return filelist

def get_target_ispcr(wildcards):
	target_dir = checkpoints.download_target.get(**wildcards).output[0]
	ISPCR_T = expand(OUTDIR + "ispcr_target/bed/{genome}.bed", 
					genome = glob_wildcards(os.path.join(target_dir, "{genome}.fna")).genome)
	return ISPCR_T

def get_offtarget_ispcr(wildcards):
	offtarget_dir = checkpoints.download_offtarget.get(**wildcards).output[0]
	ISPCR_OT = expand(OUTDIR + "ispcr_offtarget/bed/{genome}.bed", 
					genome = glob_wildcards(os.path.join(offtarget_dir, "{genome}.fna")).genome)
	return ISPCR_OT


##################################################
# Gather primer results from checkpoints
##################################################
# Get output split primers
def get_primersets_t(wildcards):
	targetdir = checkpoints.split_amplicons_target.get(**wildcards).output.alndir
	global PRIMERS_T
	PRIMERS_T = expand(OUTDIR + "ispcr_target/amplicons/{primer}.fasta",
					primer = glob_wildcards(os.path.join(targetdir, "{primer}.fasta")).primer)
	return PRIMERS_T

def get_primersets_ot(wildcards):
	offtargetdir = checkpoints.split_amplicons_offtarget.get(**wildcards).output.alndir
	global PRIMERS_OT
	PRIMERS_OT = expand(OUTDIR + "ispcr_offtarget/amplicons/{primer}.fasta",
					primer = glob_wildcards(os.path.join(offtargetdir, "{primer}.fasta")).primer)
	return PRIMERS_OT

# Get output parsed primer alignments
def get_primers_target_tsv(wildcards):
	targetdir = checkpoints.split_amplicons_target.get(**wildcards).output.alndir
	PRIMERS_T = expand(OUTDIR + "ispcr_target/primer_parsed/{primer}.tsv",
					primer = glob_wildcards(os.path.join(targetdir, "{primer}.fasta")).primer)
	return PRIMERS_T

def get_primers_offtarget_tsv(wildcards):
	offtargetdir = checkpoints.split_amplicons_offtarget.get(**wildcards).output.alndir
	PRIMERS_T = expand(OUTDIR + "ispcr_offtarget/primer_parsed/{primer}.tsv",
					primer = glob_wildcards(os.path.join(offtargetdir, "{primer}.fasta")).primer)
	return PRIMERS_T


############################################################
# Gather probe results from checkpoints
############################################################
def get_probes_t(wildcards):
	target_dir = checkpoints.parse_blast.get(**wildcards).output.probe_target_dir
	TARGET = expand(OUTDIR + "ispcr_target/probes/{probe}.fasta", 
					probe = glob_wildcards(os.path.join(target_dir, "{probe}.fasta")).probe)
	return TARGET

def get_probes_ot(wildcards):
	offtarget_dir = checkpoints.parse_blast.get(**wildcards).output.probe_offtarget_dir
	OFFTARGET = expand(OUTDIR + "ispcr_offtarget/probes/{probe}.fasta", 
					probe = glob_wildcards(os.path.join(offtarget_dir, "{probe}.fasta")).probe)
	return OFFTARGET

def get_probes_tsv_t(wildcards):
	target_dir = checkpoints.parse_blast.get(**wildcards).output.probe_target_dir
	TARGET = expand(OUTDIR + "ispcr_target/probe_parsed/{probe}.tsv", 
					probe = glob_wildcards(os.path.join(target_dir, "{probe}.fasta")).probe)
	return TARGET

def get_probes_tsv_ot(wildcards):
	offtarget_dir = checkpoints.parse_blast.get(**wildcards).output.probe_offtarget_dir
	OFFTARGET = expand(OUTDIR + "ispcr_offtarget/probe_parsed/{probe}.tsv", 
					probe = glob_wildcards(os.path.join(offtarget_dir, "{probe}.fasta")).probe)
	return OFFTARGET



if config["PROBES"] == "no":
	rule all:
		input:
			OUTDIR + "status/prepare_primers.txt",
			OUTDIR + "status/agg_target.txt",
			OUTDIR + "status/agg_offtarget.txt",
			OUTDIR + "status/collate_ispcr_target.txt",
			OUTDIR + "status/collate_ispcr_offtarget.txt",
			get_primersets_t,
			get_primersets_ot,
			OUTDIR + "status/collate_primers_target.txt",
			OUTDIR + "status/collate_primers_offtarget.txt",
			OUTDIR + "status/summary_primers_target.txt",
			OUTDIR + "status/summary_primers_offtarget.txt",
			
else:
	rule all:
		input:
			OUTDIR + "status/prepare_primers.txt",
			OUTDIR + "status/prepare_probes.txt",
			OUTDIR + "status/agg_target.txt",
			OUTDIR + "status/agg_offtarget.txt",
			OUTDIR + "status/collate_ispcr_target.txt",
			OUTDIR + "status/collate_ispcr_offtarget.txt",
			get_primersets_t,
			get_primersets_ot,
			OUTDIR + "status/collate_primers_target.txt",
			OUTDIR + "status/collate_primers_offtarget.txt",
			OUTDIR + "status/summary_primers_target.txt",
			OUTDIR + "status/summary_primers_offtarget.txt",
			get_probes_t,
			get_probes_ot,
			OUTDIR + "status/collate_probes_target.txt",
			OUTDIR + "status/collate_probes_offtarget.txt",
			OUTDIR + "status/summary_probes_target.txt",
			OUTDIR + "status/summary_probes_offtarget.txt"



include: "rules/01_get_genomes_denv.smk"
include: "rules/02_ispcr.smk"
include: "rules/03_align_primers.smk"
include: "rules/04_align_probes.smk"

# if config["USE_ASSEMBLY"] == "yes":
# 	include: "rules/03a_fastani.smk"
# else:
# 	include: "rules/03b_mafft.smk"
	
# include: "rules/04_statistics.smk"

