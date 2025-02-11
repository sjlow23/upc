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


checkpoint ispcr_target:
	input:
		genomes = GENOMES_TARGET + "{genome}.fna",
	output:
		bed = OUTDIR + "ispcr_target/bed/{genome}.bed",
		amp = OUTDIR + "ispcr_target/amplicon/{genome}.fasta",
	conda: "../envs/ispcr.yaml"
	threads: 1
	params:
		outdir = OUTDIR,
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		tile_size = TILE_SIZE,
		step_size = STEP_SIZE,
	shell:
		"""
		isPcr {input.genomes} {params.outdir}/primers.txt {output.bed} \
			-minPerfect={params.min_perfect} \
			-tileSize={params.tile_size} \
			-maxSize={params.max_size} \
			-stepSize={params.step_size} \
			-out=bed

		isPcr {input.genomes} {params.outdir}/primers.txt {output.amp} \
			-minPerfect={params.min_perfect} \
			-tileSize={params.tile_size} \
			-maxSize={params.max_size} \
			-stepSize={params.step_size} \
			-out=fa
		
		# find {params.outdir}/ispcr_target/bed -type f -name "*.bed" -empty | xargs --no-run-if-empty rm
		# find {params.outdir}/ispcr_target/amplicon -type f -name "*.fasta" -empty | xargs --no-run-if-empty rm
		
		"""


checkpoint ispcr_offtarget:
	input:
		genomes = GENOMES_OFFTARGET + "{genome}.fna",
	output:
		bed = OUTDIR + "ispcr_offtarget/bed/{genome}.bed",
		amp = OUTDIR + "ispcr_offtarget/amplicon/{genome}.fasta",
	conda: "../envs/ispcr.yaml"
	threads: 1
	params:
		outdir = OUTDIR,
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		tile_size = TILE_SIZE,
		step_size = STEP_SIZE,
	shell:
		"""
		isPcr {input.genomes} {params.outdir}/primers.txt {output.bed} \
			-minPerfect={params.min_perfect} \
			-tileSize={params.tile_size} \
			-maxSize={params.max_size} \
			-stepSize={params.step_size} \
			-out=bed

		isPcr {input.genomes} {params.outdir}/primers.txt {output.amp} \
			-minPerfect={params.min_perfect} \
			-tileSize={params.tile_size} \
			-maxSize={params.max_size} \
			-stepSize={params.step_size} \
			-out=fa

		# find {params.outdir}/ispcr_offtarget/bed -type f -name "*.bed" -empty | xargs --no-run-if-empty rm
		# find {params.outdir}/ispcr_offtarget/amplicon -type f -name "*.fasta" -empty | xargs --no-run-if-empty rm
		
		"""


rule collate_ispcr_target:
	input:
		get_target_ispcr,
	output:
		targetamp = OUTDIR + "ispcr_target/target_amplicons.fasta",
		targetdedup = OUTDIR + "ispcr_target/target_amplicons_dedup.fasta",
		status = OUTDIR + "status/collate_ispcr_target.txt",
	params:
		outdir = OUTDIR,
		targetdir = OUTDIR + "ispcr_target",
	conda: "../envs/seqkit.yaml"
	shell:
		"""
		cat {params.targetdir}/bed/*.bed > {params.targetdir}/target.bed
		cat {params.targetdir}/amplicon/*.fasta > {output.targetamp}

		sed -i "s/ /--/1" {output.targetamp}
		seqkit rmdup -s -D {params.targetdir}/duplicates_target.txt -o {output.targetdedup} {output.targetamp}
		rm -rf {params.targetdir}/bed {params.targetdir}/amplicon

		touch {output.status}
		"""


rule collate_ispcr_offtarget:
	input:
		get_offtarget_ispcr,
	output:
		offtargetamp = OUTDIR + "ispcr_offtarget/offtarget_amplicons.fasta",
		offtargetdedup = OUTDIR + "ispcr_target/offtarget_amplicons_dedup.fasta",
		status = OUTDIR + "status/collate_ispcr_offtarget.txt",
	params:
		outdir = OUTDIR,
		offtargetdir = OUTDIR + "ispcr_offtarget",
	conda: "../envs/seqkit.yaml"
	shell:
		"""
		cat {params.offtargetdir}/bed/*.bed > {params.offtargetdir}/offtarget.bed
		cat {params.offtargetdir}/amplicon/*.fasta > {output.offtargetamp}

		sed -i "s/ /--/1" {output.offtargetamp}
		seqkit rmdup -s -D {params.targetdir}/duplicates_offtarget.txt -o {output.offtargetdedup} {output.offtargetamp}
		rm -rf {params.offtargetdir}/bed {params.offtargetdir}/amplicon

		touch {output.status}
		"""