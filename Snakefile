configfile: "config/config_flu.yaml"


OUTDIR = config["OUTDIR"]
CPU = config["CPU"]
REF = config["REF"]

PRIMERS = config["PRIMERS"]
MAX_AMPLICON_SIZE = config["MAX_AMPLICON_SIZE"]
MIN_PERFECT = config["MIN_PERFECT"]
TILE_SIZE = config["TILE_SIZE"]

DOWNLOAD_TARGET = config["DOWNLOAD_TARGET"]
DOWNLOAD_OFFTARGET = config["DOWNLOAD_OFFTARGET"]
#TARGET_SP_TAXID = config["TARGET_SP_TAXID"]
#OFFTARGET_SP_TAXID = config["OFFTARGET_SP_TAXID"]
DEREP = config["DEREP"]
DATABASE = "database/"
DOMAIN = config["DOMAIN"]
ASSEMBLY_LEVEL = config["ASSEMBLY_LEVEL"]
DB = config["DB"]
USE_ASSEMBLY = config["USE_ASSEMBLY"]


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

GENOMES_T = []
GENOMES_OT = []

# Use as target in rule all to get list of target genomes
def get_target_genomes(wildcards):
	target_dir = checkpoints.download_target.get(**wildcards).output[0]
	global GENOMES_T
	GENOMES_T, = glob_wildcards(os.path.join(target_dir, "{genome}.fna"))
	filelist = expand(os.path.join(target_dir, "{genome}.fna"), genome = GENOMES_T)
	return filelist

# Use as target in rule all to get list of off-target genomes
def get_offtarget_genomes(wildcards):
	offtarget_dir = checkpoints.download_offtarget.get(**wildcards).output[0]
	global GENOMES_OT
	GENOMES_OT, = glob_wildcards(os.path.join(offtarget_dir, "{genome}.fna"))
	filelist = expand(os.path.join(offtarget_dir, "{genome}.fna"), genome = GENOMES_OT)
	return filelist

# Use as target in rule all to check primer sets 
def get_primersets_t(wildcards):
	targetdir = checkpoints.split_amplicons_target.get(**wildcards).output[0]
	global PRIMERS_T
	PRIMERS_T, = glob_wildcards(os.path.join(targetdir, "{tprimerset}.fasta"))
	primerlist_t = expand(os.path.join(targetdir, "{tprimerset}.fasta"), tprimerset = PRIMERS_T)
	return primerlist_t

def get_primersets_ot(wildcards):
	offtargetdir = checkpoints.split_amplicons_offtarget.get(**wildcards).output[0]
	global PRIMERS_OT
	PRIMERS_OT, = glob_wildcards(os.path.join(offtargetdir, "{otprimerset}.fasta"))
	primerlist_ot = expand(os.path.join(offtargetdir, "{otprimerset}.fasta"), otprimerset = PRIMERS_OT)
	return primerlist_ot

rule all:
	input:
		OUTDIR + "status/prepare_primers.txt",
		OUTDIR + "status/agg_target.txt",
		OUTDIR + "status/agg_offtarget.txt",
		OUTDIR + "status/search_target.txt",
		OUTDIR + "status/search_offtarget.txt",
		OUTDIR + "status/probesearch_target.txt",
		OUTDIR + "status/probesearch_offtarget.txt",
		OUTDIR + "status/collate.txt",
		OUTDIR + "status/split_target.txt",
		OUTDIR + "status/split_offtarget.txt",
		OUTDIR + "status/check_split.txt",
		OUTDIR + "status/align_amplicon_target.txt",
		OUTDIR + "status/align_amplicon_offtarget.txt",
		OUTDIR + "status/phylogeny.txt",
		OUTDIR + "phylogeny/metadata_subset.tsv",
		OUTDIR + "status/primerstats.txt",



include: "rules/01_get_genomes.smk"
include: "rules/02_primersearch.smk"
include: "rules/03_align.smk"

if config["USE_ASSEMBLY"] == "yes":
	include: "rules/03a_fastani.smk"
else:
	include: "rules/03b_mafft.smk"
	
include: "rules/04_statistics.smk"

