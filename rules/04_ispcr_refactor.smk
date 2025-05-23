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
		get_offtarget_genomes_tmp,
	output:
		status = OUTDIR + "status/agg_offtarget.txt"
	shell:
		"""
		touch {output.status}
		"""

checkpoint remove_overlaps:
	input:
		rules.agg_target.output.status,
		rules.agg_offtarget.output.status,
	output:
		offtargetdir = directory(GENOMES_OFFTARGET),
		status = OUTDIR + "status/remove_overlaps.txt"
	threads: 1
	params:
		outdir = OUTDIR,
		newdir = GENOMES_OFFTARGET,
		offtargetdir = GENOMES_OFFTARGET_TMP,
		target = OUTDIR + "target_genomes.txt",
		offtarget = OUTDIR + "offtarget_genomes.txt"
	shell:
		"""
		mkdir -p {params.newdir}

		if grep -wF -f {params.offtarget} {params.target}; then
			grep -wF -f {params.offtarget} {params.target} > {params.outdir}/remove
			for i in `cat {params.outdir}/remove`; do rm {params.offtargetdir}/"$i"; done
			rm {params.outdir}/remove
		fi
				
		mv {params.offtargetdir}/* {output.offtargetdir}/
		rm -rf {params.offtargetdir}

		touch {output.status}
		"""

rule reagg_offtarget:
	input:
		get_offtarget_genomes,
	output:
		status = OUTDIR + "status/reagg_offtarget.txt"
	shell:
		"""
		touch {output.status}
		"""

checkpoint ispcr_target:
	input:
		rules.prepare_primers.output.status,
		genomes = GENOMES_TARGET + "{genome}.fna",
	output:
		#bed = OUTDIR + "ispcr_target/bed/{genome}.bed",
		amp = OUTDIR + "ispcr_target/amplicon/{genome}.fasta",
	conda: "../envs/ispcr.yaml"
	threads: 2
	params:
		outdir = OUTDIR,
		primers = OUTDIR + "primers.txt",
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		tile_size = TILE_SIZE,
		step_size = STEP_SIZE,
	shell:
		"""
		isPcr {input.genomes} {params.primers} {output.amp} \
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
		rules.prepare_primers.output.status,
		genomes = GENOMES_OFFTARGET + "{genome}.fna",
	output:
		#bed = OUTDIR + "ispcr_offtarget/bed/{genome}.bed",
		amp = OUTDIR + "ispcr_offtarget/amplicon/{genome}.fasta",
	conda: "../envs/ispcr.yaml"
	threads: 2
	params:
		outdir = OUTDIR,
		primers = OUTDIR + "primers.txt",
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		tile_size = TILE_SIZE,
		step_size = STEP_SIZE,
	shell:
		"""
		isPcr {input.genomes} {params.primers} {output.amp} \
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
		status = OUTDIR + "status/collate_ispcr_target.txt",
	params:
		targetdir = OUTDIR + "ispcr_target",
		max_size = MAX_AMPLICON_SIZE
	conda: "../envs/seqkit.yaml"
	shell:
		"""
		cat {params.targetdir}/amplicon/*.fasta > {params.targetdir}/merged.fasta
		seqkit seq --max-len {params.max_size} {params.targetdir}/merged.fasta > {output.targetamp}
		rm {params.targetdir}/merged.fasta
		touch {output.status}
		"""

checkpoint missing_samples:
	input:
		rules.collate_ispcr_target.output.status,
		rules.collate_ispcr_target.output.targetamp
	output:
		outdir = directory(TARGET_NOHITS)
	params:
		outdir = OUTDIR,
		ispcrdir = OUTDIR + "ispcr_target",
		ampdir = OUTDIR + "ispcr_target/amplicons",
		resultdir = OUTDIR + "target_nohits",
		subsample = SUBSAMPLE_TARGET,
	shell:
		"""
		mkdir -p {output.outdir}
		grep ">" {params.ispcrdir}/target_amplicons.fasta | sed 's/:/\\t/1' | sed 's/>//g' | awk '{{ print $1, $3 }}' OFS="\\t" | sort | uniq > {params.resultdir}/hits.txt

		if [[ {params.subsample} == "no" ]]
		then
			genomelist="{params.outdir}/target_genomes.txt"
		else
			genomelist="{params.outdir}/target_genomes_subsampled.txt"
		fi

		for primer in `cut -f1 {params.outdir}/primers.txt | sort | uniq`; do 
			touch {params.resultdir}/"$primer".txt
			if grep -qw "$primer" {params.resultdir}/hits.txt; then
				grep -w "$primer" {params.resultdir}/hits.txt | cut -f1 | sort | uniq | awk '{{ print $0".fna" }}' >> {params.resultdir}/"$primer".txt
			fi

			#if [[ -s {params.resultdir}/"$primer".txt ]]; then 
			comm -23 "$genomelist" {params.resultdir}/"$primer".txt > {params.resultdir}/"$primer".missing
			sed -i 's/.fna//g' {params.resultdir}/"$primer".missing  
			#fi
		done	 
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
		#primers = OUTDIR + "primers.txt",
	output:
		amplicons = OUTDIR + "bbmap_amplicons/target_amplicons.fasta",
		status = OUTDIR + "status/get_missing_amplicons.txt"
	params:
		outdir = OUTDIR,
		primers = OUTDIR + "primers.txt",
		primerdir = OUTDIR + "target_nohits/",
		resultdir = OUTDIR + "bbmap_amplicons/",
		max_mismatch = MAX_MISMATCH,
		targetdir = GENOMES_TARGET,
		max_size = 200,
	conda: "../envs/align.yaml"
	threads: CPU
	shell:
		"""
		mkdir -p {params.primerdir} {params.resultdir}
		
		ls {params.primerdir}/*.txt | awk -F "/" '{{ print $NF }}' | grep -v hits | sed 's/.txt//g' > {params.outdir}/primerlist.txt			

		myprimerdir={params.primerdir}
		myresultdir={params.resultdir}
		mytargetdir={params.targetdir}

		while read primer
		do
			if [[ -s {params.primerdir}/"$primer".missing ]]
			then
				myprimer=$primer
				grep -w "$primer" {params.primers} > {params.resultdir}/input_"$primer".txt
				awk '{{ print ">"$1"\\n"$2 }}' {params.resultdir}/input_"$primer".txt > {params.primerdir}/"$primer".fwd.fasta
				awk '{{ print ">"$1"\\n"$3 }}' {params.resultdir}/input_"$primer".txt > {params.primerdir}/"$primer".rev.fasta

				cat {params.primerdir}/"$primer".missing | awk '{{ print "./scripts/extract_amplicon.sh", $0, "'$mytargetdir'", "'$myresultdir'", "'$myprimerdir'", "'$myprimer'" }}' > {params.resultdir}/run_"$primer".sh
				cat {params.resultdir}/run_"$primer".sh | parallel -j {threads}
				cat {params.resultdir}/*"$primer"_amp.fasta > {params.resultdir}/"$primer"_merged.fna
				rm {params.resultdir}/*"$primer"_amp.fasta {params.resultdir}/input_"$primer".txt {params.resultdir}/run_"$primer".sh
			fi
		done < {params.outdir}/primerlist.txt

		# Merge fasta from all primers
		cat {params.resultdir}/*_merged.fna > {params.resultdir}/merged.fasta

		# Keep only amplicons of expected size
		seqkit seq --max-len {params.max_size} {params.resultdir}/merged.fasta > {output.amplicons}

		rm {params.resultdir}/merged.fasta {params.outdir}/primerlist.txt
		rm -rf {params.primerdir}

		touch {output.status}	
		"""


rule collate_ispcr_bbmap:
	input:
		ispcr = rules.collate_ispcr_target.output.targetamp,
		bbmap = rules.get_missing_amplicons.output.amplicons,
	output:
		merged = OUTDIR + "ispcr_target/merged_amplicons.fasta",
		targetdedup = OUTDIR + "ispcr_target/merged_amplicons_dedup.fasta",
		status = OUTDIR + "status/collate_ispcr_bbmap.txt"
	params:
		outdir = OUTDIR,
		targetdir = OUTDIR + "ispcr_target"
	conda: "../envs/seqkit.yaml"
	shell:
		"""
		cat {input.ispcr} {input.bbmap} | seqkit seq --upper-case > {output.merged}
		cp {output.merged} {output.targetdedup}
		
		rm -rf {params.targetdir}/amplicon
		touch {output.status}
		"""

rule collate_ispcr_offtarget:
	input:
		get_offtarget_ispcr,
	output:
		offtargetamp = OUTDIR + "ispcr_offtarget/offtarget_amplicons.fasta",
		offtargetdedup = OUTDIR + "ispcr_offtarget/offtarget_amplicons_dedup.fasta",
		status = OUTDIR + "status/collate_ispcr_offtarget.txt",
	params:
		outdir = OUTDIR,
		offtargetdir = OUTDIR + "ispcr_offtarget",
	conda: "../envs/seqkit.yaml"
	shell:
		"""
		#cat {params.offtargetdir}/bed/*.bed > {params.offtargetdir}/offtarget.bed
		cat {params.offtargetdir}/amplicon/*.fasta | seqkit seq --upper-case > {output.offtargetamp}

		#sed -i 's/ /--/1' {output.offtargetamp}
		seqkit rmdup -s -D {params.offtargetdir}/duplicates_offtarget.txt -o {output.offtargetdedup} {output.offtargetamp}
		
		if [[ -s {params.offtargetdir}/duplicates_offtarget.txt ]]; then
			cut -f2 {params.offtargetdir}/duplicates_offtarget.txt | sed 's/, /\\t/1' | \
				awk -F "\\t" '{{ print $2, $1, $1 }}' OFS="\t" | \
				sed 's/\\t/, /1' > {params.offtargetdir}/duplicates_offtarget.tab
		fi
		rm -rf {params.offtargetdir}/bed {params.offtargetdir}/amplicon

		if [[ -s {output.offtargetdedup} ]]
		then
			echo "Present" > {output.status}
		else
			touch {output.status}
		fi
		"""


