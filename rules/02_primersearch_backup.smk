rule agg_target:
	input:
		get_target_genomes,
	output:
		status = OUTDIR + "status/agg_target.txt",
	shell:
		"""
		touch {output.status}
		"""

rule agg_offtarget:
	input:
		get_offtarget_genomes,
	output:
		status = OUTDIR + "status/agg_offtarget.txt"
	shell:
		"""
		touch {output.status}
		"""


rule search_target:
	input:
		status = rules.agg_target.output.status,
		genomes = GENOMES_TARGET + "{genome}.fna",
		
	output:
		bed = OUTDIR + "ispcr_target/output_bed/{genome}.bed",
		amplicons = OUTDIR + "ispcr_target/output_fasta/{genome}.fasta",
	
	conda: "../envs/ispcr.yaml"

	threads: 8

	params:
		outdir = OUTDIR,
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		genome = "{genome}",
		seed = TILE_SIZE,

	shell:
		"""
		isPcr {input.genomes} {params.outdir}/primers.txt {output.bed} \
			-minPerfect={params.min_perfect} \
			-tileSize={params.seed} \
			-maxSize={params.max_size} \
			-out=bed

		isPcr {input.genomes} {params.outdir}/primers.txt {output.amplicons} \
			-minPerfect={params.min_perfect} \
			-tileSize={params.seed} \
			-maxSize={params.max_size} \
			-out=fa
		
		"""

def get_nonempty_files(wildcards, input):
    return [i for i in input if os.path.getsize(i)>0]

def get_target_genomes(wildcards):
	target_dir = checkpoints.download_target.get(**wildcards).output[0]
	#global GENOMES_T
	GENOMES_T, = glob_wildcards(os.path.join(target_dir, "{genome}.fna"))
	#filelist = expand(os.path.join(target_dir, "{genome}.fna"), genome = GENOMES_T)
	return GENOMES_T

get_target_genomes

rule get_nonempty_target_amplicon:
	input:
		expand(os.path.join(OUTDIR, "ispcr_target/output_fasta/{genome}.fasta"),	genome = GENOMES_T)
	output:
		ampdir = OUTDIR + "ispcr_target/amplicons/",
	params:
		nonempty = get_nonempty_files
	shell:
		"""
		mv {params.nonempty} {output.ampdir}
		"""

rule get_nonempty_target_bed:
	input:
		get_target_genomes,
		expand(os.path.join(OUTDIR, "ispcr_target/output_bed/{genome}.bed"), genome = GENOMES_T)
	output:
		beddir = OUTDIR + "ispcr_target/bed/",
	params:
		nonempty = get_nonempty_files
	shell:
		"""
		mv {params.nonempty} {output.beddir}
		"""


rule search_offtarget:
	input:
		rules.agg_offtarget.output.status,
		genomes = GENOMES_OFFTARGET + "{genome}.fna",
			
	output:
		bed = OUTDIR + "ispcr_offtarget/output_bed/{genome}.bed",
		amplicons = OUTDIR + "ispcr_offtarget/output_fasta/{genome}.fasta",
	
	conda: "../envs/ispcr.yaml"

	threads: 8

	params:
		outdir = OUTDIR,
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		genome = "{genome}",
		seed = TILE_SIZE,
	
	shell:
		"""
		isPcr {input.genomes} {params.outdir}/primers.txt {output.bed} \
			-minPerfect={params.min_perfect} \
			-tileSize={params.seed} \
			-maxSize={params.max_size} \
			-out=bed

		isPcr {input.genomes} {params.outdir}/primers.txt {output.amplicons} \
			-minPerfect={params.min_perfect} \
			-tileSize={params.seed} \
			-maxSize={params.max_size} \
			-out=fa
		
		"""

rule probesearch_target_amplicon:
	input:
		amplicons = OUTDIR + "ispcr_target/amplicons/{genome}.fasta",

	output:
		sam_perfect = OUTDIR + "ispcr_target/probes/{genome}_nomm.fasta",
		sam_mm = OUTDIR + "ispcr_target/probes/{genome}_2mm.fasta",

	conda: "../envs/ispcr.yaml"

	threads: 4

	params:
		genome = "{genome}",
		outdir = OUTDIR

	shell:
		"""
		bbmap.sh \
		in={params.outdir}/probes.fasta \
		ref={params.nonempty} \
		nodisk \
		noheader=t \
		ambig=all \
		vslow \
		perfectmode \
		maxsites=100000 \
		threads={threads} \
		outm={output.sam_perfect}

		bbmap.sh \
		in={params.outdir}/probes.fasta \
		ref={input.amplicons} \
		nodisk \
		noheader=t \
		ambig=all \
		vslow \
		editfilter=1 \
		maxsites=100000 \
		threads={threads} \
		outm={output.sam_mm}

		"""

#Blast against per-genome amplicons, need to match primer name to probe name
#Blast against whole genome

rule collate_results:
	input:
		get_target_output,
		get_offtarget_output,
		get_probe_target_output

	output:
		bedfile_target = OUTDIR + "ispcr_target/target.bed",
		bedfile_offtarget = OUTDIR + "ispcr_offtarget/offtarget.bed",
		fasta_target = OUTDIR + "ispcr_target/target_amplicons.fasta",
		fasta_offtarget = OUTDIR + "ispcr_offtarget/offtarget_amplicons.fasta",
		status = OUTDIR + "status/collate.txt"

	params:
		target_bed = OUTDIR + "ispcr_target/output_bed",
		offtarget_bed = OUTDIR + "ispcr_offtarget/output_bed",
		target_fasta = OUTDIR + "ispcr_target/output_fasta",
		offtarget_fasta = OUTDIR + "ispcr_offtarget/output_fasta",
		targetdir = OUTDIR + "ispcr_target",
		offtargetdir = OUTDIR + "ispcr_offtarget",

	shell:
		"""
		cat {params.target_bed}/*.bed > {output.bedfile_target}
		cat {params.offtarget_bed}/*.bed > {output.bedfile_offtarget}
		
		cat {params.target_fasta}/*.fasta > {output.fasta_target}
		cat {params.offtarget_fasta}/*.fasta > {output.fasta_offtarget}

		rm -rf {params.target_bed} {params.offtarget_bed}
		rm -rf {params.target_fasta} {params.offtarget_fasta}

		touch {output.status}

		"""