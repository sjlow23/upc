checkpoint missing_samples:
	input:
		rules.collate_ispcr_target.output.status,
		rules.collate_ispcr_target.output.targetamp
	output:
		outdir = directory(TARGET_NOHITS)
	params:
		outdir = OUTDIR,
		ampdir = OUTDIR + "ispcr_target/amplicons/",
		resultdir = OUTDIR + "target_nohits/",
		subsample = SUBSAMPLE_TARGET,
	shell:
		"""
		grep ">" target_amplicons.fasta | sed 's/:/\t/1' | sed 's/>//g' | awk '{{ print $1, $3 }}' OFS="\\t" > {params.resultdir}/hits.txt

		mkdir -p {output.outdir}

		if [[ {params.subsample}  != "no" ]]
		then
			genomelist="{params.outdir}/target_genomes_subsampled.txt"
		else
			genomelist="{params.outdir}/target_genomes.txt"
		fi

		for primer in `cut -f1 {params.resultdir}/primers.txt | sort | uniq`; do \
			grep -w "$primer" {params.resultdir}/hits.txt | cut -f1 | sort | uniq > {params.resultdir}/"$primer".txt
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
		amplicons = OUTDIR + "bbmap_amplicons/target_amplicons.fasta",
		status = OUTDIR + "status/get_missing_amplicons.txt"
	params:
		outdir = OUTDIR,
		primerdir = OUTDIR + "target_nohits/",
		resultdir = OUTDIR + "bbmap_amplicons/",
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
				noheader=t \
				trimreaddescriptions=t \
				out={params.resultdir}/"$genome".fwd.sam \
				ref={params.primerdir}/"$primer".fwd cutoff=0.7;

				msa.sh in={params.targetdir}/"$genome".fna \
				noheader=t \
				trimreaddescriptions=t \
				out={params.resultdir}/"$genome".rev.sam \
				rcomp=t \
				ref={params.primerdir}/"$primer".rev cutoff=0.7

				cutprimers.sh \
				trimreaddescriptions=t \
				sam1={params.resultdir}/"$genome".fwd.sam \
				sam2={params.resultdir}/"$genome".rev.sam \
				in={params.targetdir}/"$genome".fna \
				out={params.resultdir}/"$genome".fasta \
				include=t

				cat {params.resultdir}/"$genome".fasta >> {params.resultdir}/"$primer"_merged.fasta

				# Add alignment info to lookup file for renaming header
				awk -F "\\t" '{{ print $3, $1, $10 }}' OFS="\t" {params.resultdir}/"$genome".fwd.sam | sed 's/\\t/--/1' >> {params.resultdir}/fwd
				awk -F "\\t" '{{ print $3, $1, $10 }}' OFS="\t" {params.resultdir}/"$genome".rev.sam | sed 's/\\t/--/1' >> {params.resultdir}/rev
				
				rm {params.resultdir}/"$genome".fwd.sam {params.resultdir}/"$genome".rev.sam {params.resultdir}/"$genome".fasta
			done

				# Rename headers in fasta
				sort -k1,1 {params.resultdir}/fwd > {params.resultdir}/fwd.sorted
				sort -k1,1 {params.resultdir}/rev > {params.resultdir}/rev.sorted
				join -1 1 -2 1 -a 1 -a 2 -e NA -o 1.1,1.2,2.2 {params.resultdir}/fwd.sorted {params.resultdir}/rev.sorted | \
				sed 's/--/\\t/1' | \
				awk '{{ print $1, $1":0+0", $2, "X", $3, $4 }}' OFS="\\t" | \
				sed 's/\\t/ /g' | \
				sed 's/ /\\t/1' > {params.resultdir}/headers.txt
				rm {params.resultdir}/fwd {params.resultdir}/rev {params.resultdir}/fwd.sorted {params.resultdir}/rev.sorted

			seqkit replace --kv-file {params.resultdir}/headers.txt --pattern "^(\\S+)" --replacement "{{kv}}" --out-file {params.resultdir}/"$primer"_amp.fasta {params.resultdir}/"$primer"_merged.fasta

			rm {params.primerdir}/"$primer".fwd {params.primerdir}/"$primer".rev
		done < {params.outdir}/primerlist.txt

		# Merge fasta from all primers
		cat {params.resultdir}/*_amp.fasta > {params.resultdir}/merged.fasta

		# Keep only amplicons of expected size
		seqkit seq --max-len {params.max_size} {params.resultdir}/merged.fasta > {output.amplicons}

		rm {params.resultdir}/*_amp.fasta {params.resultdir}/merged.fasta
		touch {output.status}
		"""


