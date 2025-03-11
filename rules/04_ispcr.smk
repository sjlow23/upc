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

		grep -wF -f {params.offtarget} {params.target} > {params.outdir}/remove

		if [[ -s {params.outdir}/remove ]]
		then
			for i in `cat {params.outdir}/remove`; do rm {params.offtargetdir}/"$i"; done
		fi
				
		mv {params.offtargetdir}/* {output.offtargetdir}/
		rm {params.outdir}/remove
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
		genomes = GENOMES_TARGET + "{genome}.fna",
	output:
		#bed = OUTDIR + "ispcr_target/bed/{genome}.bed",
		amp = OUTDIR + "ispcr_target/amplicon/{genome}.fasta",
	conda: "../envs/ispcr.yaml"
	threads: 2
	params:
		outdir = OUTDIR,
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		tile_size = TILE_SIZE,
		step_size = STEP_SIZE,
	shell:
		"""
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
		#bed = OUTDIR + "ispcr_offtarget/bed/{genome}.bed",
		amp = OUTDIR + "ispcr_offtarget/amplicon/{genome}.fasta",
	conda: "../envs/ispcr.yaml"
	threads: 2
	params:
		outdir = OUTDIR,
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		tile_size = TILE_SIZE,
		step_size = STEP_SIZE,
	shell:
		"""
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

		if [[ {params.subsample} == "" || {params.subsample} == "no" ]]
		then
			genomelist="{params.outdir}/target_genomes.txt"
		else
			genomelist="{params.outdir}/target_genomes_subsampled.txt"
		fi

		for primer in `cut -f1 {params.outdir}/primers.txt | sort | uniq`; do \
			grep -w "$primer" {params.resultdir}/hits.txt | cut -f1 | sort | uniq | awk '{{ print $0".fna" }}' > {params.resultdir}/"$primer".txt
			comm -23 $genomelist {params.resultdir}/"$primer".txt > {params.resultdir}/"$primer".missing; done
	
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
		max_size = 200,
	conda: "../envs/align.yaml"
	threads: 4
	shell:
		"""
		mkdir -p {params.primerdir} {params.resultdir}
		
		ls {params.primerdir}/*.txt | awk -F "/" '{{ print $NF }}' | grep -v hits | sed 's/.txt//g' > {params.outdir}/primerlist.txt

		while read primer
		do
			grep -w "$primer" {input.primers} > input.txt
			awk '{{ print ">"$1"\\n"$2 }}' input.txt > {params.primerdir}/"$primer".fwd.fasta
			awk '{{ print ">"$1"\\n"$3 }}' input.txt > {params.primerdir}/"$primer".rev.fasta

			for genome in `cat {params.primerdir}/"$primer".missing | sed 's/.fna//g'`; do \
				msa.sh in={params.targetdir}/"$genome".fna \
				noheader=t \
				trimreaddescriptions=t \
				out={params.resultdir}/"$genome".fwd.sam \
				ref={params.primerdir}/"$primer".fwd.fasta cutoff=0.75

				msa.sh in={params.targetdir}/"$genome".fna \
				noheader=t \
				trimreaddescriptions=t \
				addr=t \
				out={params.resultdir}/"$genome".rev.sam \
				ref={params.primerdir}/"$primer".rev.fasta cutoff=0.75

				cutprimers.sh \
				trimreaddescriptions=t \
				sam1={params.resultdir}/"$genome".fwd.sam \
				sam2={params.resultdir}/"$genome".rev.sam \
				in={params.targetdir}/"$genome".fna \
				out={params.resultdir}/"$genome".fasta \
				include=t \
				fake=f

				for i in {params.resultdir}/*.sam; do \
				awk -F "\t" 'BEGIN {{
					comp["A"] = "T"; comp["T"] = "A"; comp["C"] = "G"; comp["G"] = "C"
				}}
				$2 != 4 && $2 != 20
				{{
					if ($1 ~ /^r_/) {{
						seq = $10
						rev_comp = ""
						for (i = length(seq); i > 0; i--) {{
							rev_comp = rev_comp comp[substr(seq, i, 1)]
						}}
						$10 = rev_comp
					}}
					print
				}}' OFS="\t" "$i" | sed 's/^r_//g' > "$i".rc; done

				cat {params.resultdir}/"$genome".fasta >> {params.resultdir}/"$primer"_merged.fasta

				# Add alignment info to lookup file for renaming header
				awk -F "\\t" '{{ print $3, $1, $4, $10 }}' OFS="\\t" {params.resultdir}/"$genome".fwd.sam.rc | sed 's/\\t/--/1' >> {params.resultdir}/fwd
				awk -F "\\t" '{{ print $3, $1, $4+length($10)-1, $10 }}' OFS="\\t" {params.resultdir}/"$genome".rev.sam.rc | sed 's/\\t/--/1' >> {params.resultdir}/rev
				
				#rm {params.resultdir}/"$genome".fwd.sam* {params.resultdir}/"$genome".rev.sam* {params.resultdir}/"$genome".fasta
			done

				# Rename headers in fasta
				sort -k1,1 {params.resultdir}/fwd > {params.resultdir}/fwd.sorted
				sort -k1,1 {params.resultdir}/rev > {params.resultdir}/rev.sorted

				join -1 1 -2 1 -a 1 -a 2 -e NA -o 1.1,1.2,2.2,1.3,2.3 {params.resultdir}/fwd.sorted {params.resultdir}/rev.sorted | \
				sed 's/--/\\t/1' | \
				awk '{{ print $1, $1":"$3"+"$4, $2, "X", $5, $6 }}' OFS="\\t" | \
				sed 's/\\t/ /g' | \
				sed 's/ /\\t/1' > {params.resultdir}/headers.txt
				#rm {params.resultdir}/fwd {params.resultdir}/rev {params.resultdir}/fwd.sorted {params.resultdir}/rev.sorted

			seqkit replace --kv-file {params.resultdir}/headers.txt --pattern "^(\\S+)" --replacement "{{kv}}" --out-file {params.resultdir}/"$primer"_amp.fasta {params.resultdir}/"$primer"_merged.fasta

			#rm {params.primerdir}/"$primer".fwd.fasta {params.primerdir}/"$primer".rev.fasta {params.resultdir}/headers.txt
		done < {params.outdir}/primerlist.txt

		# Merge fasta from all primers
		cat {params.resultdir}/*amp.fasta > {params.resultdir}/merged.fasta

		# Keep only amplicons of expected size
		seqkit seq --max-len {params.max_size} {params.resultdir}/merged.fasta > {output.amplicons}

		#rm {params.resultdir}/*_amp.fasta {params.resultdir}/merged.fasta {params.outdir}/primerlist.txt 
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
		cat {input.ispcr} {input.bbmap} > {output.merged}

		seqkit rmdup -s -D {params.targetdir}/duplicates_target.txt -o {output.targetdedup} {output.merged}
		
		if [[ -s {params.targetdir}/duplicates_target.txt ]]; then
			cut -f2 {params.targetdir}/duplicates_target.txt | sed 's/, /\\t/1' | \
				awk -F "\\t" '{{ print $2, $1, $1 }}' OFS="\\t" | \
				sed 's/\\t/, /1' > {params.targetdir}/duplicates_target.tab
		fi
		
		rm -rf {params.targetdir}/bed {params.targetdir}/amplicon
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
		cat {params.offtargetdir}/amplicon/*.fasta > {output.offtargetamp}

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


