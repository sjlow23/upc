checkpoint missing_samples:
	input:
		get_primersets_t,
	output:
		outdir = directory(TARGET_NOHITS),
	params:
		outdir = OUTDIR,
		ampdir = OUTDIR + "ispcr_target/amplicons/",
		resultdir = OUTDIR + "target_nohits/",
		subsample = SUBSAMPLE_TARGET,
	shell:
		"""
		mkdir -p {output.outdir}

		if [[ {params.subsample}  != "no" ]]
		then
			genomelist="{params.outdir}/target_genomes_subsampled.txt"
		else
			genomelist="{params.outdir}/target_genomes.txt"
		fi

		for i in `ls {params.ampdir}/*.fasta | awk -F "/" '{{ print $NF }}' | sed 's/.fasta//g'`; do \
			grep ">" {params.ampdir}/"$i".fasta | awk -F ":" '{{ print $1 }}' | sed 's/>//g' | awk '{{ print $0".fna" }}' | sort | uniq > {params.resultdir}/"$i".txt
			comm -23 $genomelist {params.resultdir}/"$i".txt > {params.resultdir}/"$i".missing; done
		"""


rule agg_missing_samples:
	input:
		get_missing_samples
	output:
		status = OUTDIR + "status/agg_missing_samples.txt"
	shell:
		"""
		touch {output.status}
		"""


rule get_missing_amplicons:
	input:
		rules.agg_missing_samples.output.status,
		primers = OUTDIR + "primers.txt",
	output:
		amplicons = OUTDIR + "seqkit_amplicons/collated_amplicons.fasta",
		status = OUTDIR + "status/get_missing_amplicons.txt"
	params:
		outdir = OUTDIR,
		primerdir = OUTDIR + "target_nohits/",
		resultdir = OUTDIR + "seqkit_amplicons/",
		max_mismatch = MAX_MISMATCH,
		targetdir = GENOMES_TARGET,
		max_size = MAX_AMPLICON_SIZE,
	conda: "../envs/align.yaml"
	threads: 4
	shell:
		"""
		mkdir -p {params.primerdir}

		ls {params.primerdir}/*.txt | awk -F "/" '{{ print $NF }}' | sed 's/.txt//g' > {params.outdir}/primerlist.txt

		while read primer
		do
			grep -w "$primer" {input.primers} > input.txt
			awk '{{ print ">"$1"\\n"$2 }}' input.txt > {params.primerdir}/"$primer".fwd
			awk '{{ print ">"$1"\\n"$3 }}' input.txt > {params.primerdir}/"$primer".rev

			for genome in `cat {params.primerdir}/"$primer".missing | sed 's/.fna//g'`; do \
				msa.sh in={params.targetdir}/"$genome".fna \
				out={params.resultdir}/"$genome".fwd.sam \
				ref={params.primerdir}/"$primer".fwd cutoff=0.7;
				msa.sh in={params.targetdir}/"$genome".fna \
				out={params.resultdir}/"$genome".rev.sam \
				rcomp=t \
				ref={params.primerdir}/"$primer".rev cutoff=0.7
			
				cutprimers.sh \
				sam1={params.resultdir}/"$genome".fwd.sam \
				sam2={params.resultdir}/"$genome".rev.sam \
				in={params.targetdir}/"$genome".fna \
				out={params.resultdir}/"$genome".fasta \
				include=t

				cat {params.resultdir}/"$genome".fasta >> {params.resultdir}/"$primer"_amp.fasta
				rm {params.resultdir}/"$genome".fwd.sam {params.resultdir}/"$genome".rev.sam {params.resultdir}/"$genome".fasta
			done
			
			rm {params.primerdir}/"$primer".fwd {params.primerdir}/"$primer".rev
		done < {params.outdir}/primerlist.txt

		cat {params.resultdir}/*_amp.fasta > {params.resultdir}/merged.fasta

		# Keep only amplicons of expected size
		seqkit seq --max-len {params.max_size} {params.resultdir}/merged.fasta > {output.amplicons}

		rm {params.resultdir}/*_amp.fasta
		touch {output.status}
		"""


		