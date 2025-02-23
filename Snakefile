configfile: "config/config_denv.yaml"

OUTDIR = directory(config["params"]["OUTDIR"])
CPU = config["params"]["CPU"]
SUBSAMPLE_TARGET = config["params"]["SUBSAMPLE_TARGET"]
SUBSAMPLE_OFFTARGET = config["params"]["SUBSAMPLE_OFFTARGET"]
PRIMERS = config["params"]["PRIMERS"]

MAX_AMPLICON_SIZE = config["ispcr"]["MAX_AMPLICON_SIZE"]
MIN_PERFECT = config["ispcr"]["MIN_PERFECT"]
TILE_SIZE = config["ispcr"]["TILE_SIZE"]
STEP_SIZE = config["ispcr"]["STEP_SIZE"]
MIN_GOOD = config["ispcr"]["MIN_GOOD"]

DOMAIN = config["datasets"]["DOMAIN"]
ASSEMBLY_LEVEL = config["datasets"]["ASSEMBLY_LEVEL"]
DB = config["datasets"]["DB"]
USE_ASSEMBLY = config["datasets"]["USE_ASSEMBLY"]
IS_FLU = config["datasets"]["IS_FLU"]

if "PROBES" in config.get("params", {}):
	PROBES = config["params"]["PROBES"]
else:
	# Check if the number of columns in PRIMERS is 4
	with open(PRIMERS, "r") as f:
		first_line = f.readline().strip()
		num_columns = len(first_line.split("\t"))  
	if num_columns == 4:
		PROBES = "yes"
	else:
		PROBES = "no"

DOWNLOAD_TARGET = config["download"]["DOWNLOAD_TARGET"]
DOWNLOAD_OFFTARGET = config["download"]["DOWNLOAD_OFFTARGET"]

if "USER_TARGET" in config.get("custom", {}):
    USER_TARGET = config["custom"]["USER_TARGET"]
if "USER_OFFTARGET" in config.get("custom", {}):
	USER_OFFTARGET = config["custom"]["USER_OFFTARGET"]
if "TARGET_SP_TAXID" in config.get("datasets", {}):
	TARGET_SP_TAXID = config["datasets"]["TARGET_SP_TAXID"]
if "OFFTARGET_SP_TAXID" in config.get("datasets", {}):
	OFFTARGET_SP_TAXID = config["datasets"]["OFFTARGET_SP_TAXID"]
# if "TARGET_TAXID" in config:
# 	TARGET_TAXID = config["TARGET_TAXID"]
# if "OFFTARGET_TAXID" in config:
# 	OFFTARGET_TAXID = config["OFFTARGET_TAXID"]



GENOMES_TARGET = OUTDIR + "genomes_target/"
GENOMES_OFFTARGET = OUTDIR + "genomes_offtarget/"
GENOMES_OFFTARGET_TMP = OUTDIR + "genomes_offtarget_tmp/"

# Use as target in rule all to get list of target genomes
def get_target_genomes(wildcards):
	target_dir = checkpoints.download_target.get(**wildcards).output[0]
	filelist = expand(os.path.join(target_dir, "{genome}.fna"), 
					genome = glob_wildcards(os.path.join(target_dir, "{genome}.fna")).genome)
	return filelist

# Use as target in rule all to get list of off-target genomes
def get_offtarget_genomes_tmp(wildcards):
	offtarget_dir = checkpoints.download_offtarget.get(**wildcards).output[0]
	filelist = expand(os.path.join(offtarget_dir, "{genome}.fna"), 
					genome = glob_wildcards(os.path.join(offtarget_dir, "{genome}.fna")).genome)
	return filelist

# Get offtarget genomes after removal of overlaps
def get_offtarget_genomes(wildcards):
	offtarget_dir = checkpoints.remove_overlaps.get(**wildcards).output[0]
	filelist = expand(os.path.join(offtarget_dir, "{genome}.fna"), 
					genome = glob_wildcards(os.path.join(offtarget_dir, "{genome}.fna")).genome)
	return filelist


def get_target_ispcr(wildcards):
	target_dir = checkpoints.download_target.get(**wildcards).output[0]
	ISPCR_T = expand(OUTDIR + "ispcr_target/amplicon/{genome}.fasta", 
					genome = glob_wildcards(os.path.join(target_dir, "{genome}.fna")).genome)
	return ISPCR_T

def get_offtarget_ispcr(wildcards):
	offtarget_dir = checkpoints.remove_overlaps.get(**wildcards).output[0]
	ISPCR_OT = expand(OUTDIR + "ispcr_offtarget/amplicon/{genome}.fasta", 
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
	target_dir = checkpoints.parse_blast_target.get(**wildcards).output.probe_target_dir
	TARGET = expand(OUTDIR + "ispcr_target/probes/{probe}.fasta", 
					probe = glob_wildcards(os.path.join(target_dir, "{probe}.fasta")).probe)
	return TARGET

def get_probes_ot(wildcards):
	offtarget_dir = checkpoints.parse_blast_offtarget.get(**wildcards).output.probe_offtarget_dir
	OFFTARGET = expand(OUTDIR + "ispcr_offtarget/probes/{probe}.fasta", 
					probe = glob_wildcards(os.path.join(offtarget_dir, "{probe}.fasta")).probe)
	return OFFTARGET

def get_probes_tsv_t(wildcards):
	target_dir = checkpoints.parse_blast_target.get(**wildcards).output.probe_target_dir
	TARGET = expand(OUTDIR + "ispcr_target/probe_parsed/{probe}.tsv", 
					probe = glob_wildcards(os.path.join(target_dir, "{probe}.fasta")).probe)
	return TARGET

def get_probes_tsv_ot(wildcards):
	offtarget_dir = checkpoints.parse_blast_offtarget.get(**wildcards).output.probe_offtarget_dir
	OFFTARGET = expand(OUTDIR + "ispcr_offtarget/probe_parsed/{probe}.tsv", 
					probe = glob_wildcards(os.path.join(offtarget_dir, "{probe}.fasta")).probe)
	return OFFTARGET



# Define common files that are always included in rule_all
rule_all_input = [
	OUTDIR + "status/prepare_primers.txt",
	OUTDIR + "status/agg_target.txt",
	OUTDIR + "status/collate_ispcr_target.txt",
	get_primersets_t,
	OUTDIR + "status/collate_primers_target.txt",
	OUTDIR + "status/summary_primers_target.txt",
]

# Handle "PROBES" condition to include probe-related files only if PROBES != "no"
if config["params"]["PROBES"] != "no" and not config["datasets"]["OFFTARGET_SP_TAXID"]:
	rule_all_input.extend([
		prepare_probes,
		OUTDIR + "status/prepare_probes.txt",
		get_probes_t,
		OUTDIR + "status/collate_probes_target.txt",
		OUTDIR + "status/summary_probes_target.txt"
	
	])


# Handle "GENOMES_OFFTARGET" condition to include primer-related files only if config["OFFTARGET_SP_TAXID"]
if config["datasets"]["OFFTARGET_SP_TAXID"]:
	rule_all_input.extend([
	OUTDIR + "status/agg_offtarget.txt",
	OUTDIR + "status/reagg_offtarget.txt",
	OUTDIR + "status/collate_ispcr_offtarget.txt",
	get_primersets_ot,
	OUTDIR + "status/collate_primers_offtarget.txt",
	OUTDIR + "status/summary_primers_offtarget.txt"
	])
	if config["params"]["PROBES"] != "no":
		rule_all_input.extend([
			get_probes_ot,
			OUTDIR + "status/collate_probes_target.txt",
			OUTDIR + "status/summary_probes_target.txt"
		])

rule all:
	input:
		rule_all_input





include: "rules/01_get_genomes_target.smk"
include: "rules/02_get_genomes_offtarget.smk"
include: "rules/03_prepare_primers_probes.smk"
include: "rules/04_ispcr.smk"
include: "rules/05_align_primers.smk"
include: "rules/06_align_probes.smk"
