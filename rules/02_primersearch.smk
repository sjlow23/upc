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
		genomes = OUTDIR + "genomes_target.txt",
	output:
		beddir = directory(OUTDIR + "ispcr_target/bed"),
		ampdir = directory(OUTDIR + "ispcr_target/amplicon"),
		status = OUTDIR + "status/search_target.txt",
	conda: "../envs/ispcr.yaml"
	threads: 8
	params:
		outdir = OUTDIR,
		genomedir = GENOMES_TARGET,
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		seed = TILE_SIZE,
	shell:
		"""
		mkdir -p {output.beddir}
		mkdir -p {output.ampdir}

		cat {input.genomes} | \
			parallel -j {threads} \
			'isPcr {params.genomedir}/"{{}}".fna {params.outdir}/primers.txt {output.beddir}/"{{}}".bed \
			-minPerfect={params.min_perfect} \
			-tileSize={params.seed} \
			-maxSize={params.max_size} \
			-out=bed'
		
		cat {input.genomes} | \
			parallel -j {threads} \
			'isPcr {params.genomedir}/"{{}}".fna {params.outdir}/primers.txt {output.ampdir}/"{{}}".fasta \
			-minPerfect={params.min_perfect} \
			-tileSize={params.seed} \
			-maxSize={params.max_size} \
			-out=fa'
		
		find {output.beddir} -type f -name "*.bed" -empty | xargs --no-run-if-empty rm
		find {output.ampdir} -type f -name "*.fasta" -empty | xargs --no-run-if-empty rm

		if [[ "$(ls -A {output.beddir})" ]]
		then
			find {output.beddir}/*.bed | awk -F "/" '{{ print $NF }}' | sed 's/.bed//g' > {params.outdir}/genomes_target_primerhits.txt
		fi

		touch {output.status}
		"""



rule search_offtarget:
	input:
		status = rules.agg_offtarget.output.status,
		genomes = OUTDIR + "genomes_offtarget.txt",
	output:
		beddir = directory(OUTDIR + "ispcr_offtarget/bed"),
		ampdir = directory(OUTDIR + "ispcr_offtarget/amplicon"),
		status = OUTDIR + "status/search_offtarget.txt",
	conda: "../envs/ispcr.yaml"
	threads: 8
	params:
		outdir = OUTDIR,
		genomedir = GENOMES_OFFTARGET,
		min_perfect = MIN_PERFECT,
		max_size = MAX_AMPLICON_SIZE,
		seed = TILE_SIZE,
	shell:
		"""
		mkdir -p {output.beddir}
		mkdir -p {output.ampdir}

		cat {input.genomes} | \
			parallel -j {threads} \
			'isPcr {params.genomedir}/"{{}}".fna {params.outdir}/primers.txt {output.beddir}/"{{}}".bed \
			-minPerfect={params.min_perfect} \
			-tileSize={params.seed} \
			-maxSize={params.max_size} \
			-out=bed'
		
		cat {input.genomes} | \
			parallel -j {threads} \
			'isPcr {params.genomedir}/"{{}}".fna {params.outdir}/primers.txt {output.ampdir}/"{{}}".fasta \
			-minPerfect={params.min_perfect} \
			-tileSize={params.seed} \
			-maxSize={params.max_size} \
			-out=fa'
		
		find {output.beddir} -type f -name "*.bed" -empty | xargs --no-run-if-empty rm
		find {output.ampdir} -type f -name "*.fasta" -empty | xargs --no-run-if-empty rm

		if [[ "$(ls -A {output.beddir})" ]]
		then
			find {output.beddir}/*.bed | awk -F "/" '{{ print $NF }}' | sed 's/.bed//g' > {params.outdir}/genomes_offtarget_primerhits.txt
		fi

		touch {output.status}
		"""


rule probesearch_amplicon_target:
	input:
		target_genomes = rules.search_target.output.status,
	output:
		status = OUTDIR + "status/probesearch_target.txt",
	conda: "../envs/ispcr.yaml"
	threads: 4
	params:
		ampdir = OUTDIR + "ispcr_target/amplicon",
		probedir_nomm = directory(OUTDIR + "ispcr_target/probes/no_mismatch"),
		probedir_2mm = directory(OUTDIR + "ispcr_target/probes/2_mismatch"),
		ispcrdir = directory(OUTDIR + "ispcr_target"),
		outdir = OUTDIR
	shell:
		"""
		if [[ -s {params.outdir}/genomes_target_primerhits.txt ]]
		then
			cat {params.outdir}/genomes_target_primerhits.txt | \
			parallel -j {threads} \
			'bbmap.sh \
			in={params.outdir}/probes.fasta \
			ref={params.ampdir}/"{{}}".fasta \
			nodisk \
			noheader=t \
			ambig=all \
			vslow \
			perfectmode \
			maxsites=100000 \
			overwrite=f \
			threads={threads} \
			outm={params.probedir_nomm}/"{{}}"_nomm.sam'

			cat {params.outdir}/genomes_target_primerhits.txt | \
			parallel -j {threads} \
			'bbmap.sh \
			in={params.outdir}/probes.fasta \
			ref={params.ampdir}/"{{}}".fasta \
			nodisk \
			noheader=t \
			ambig=all \
			vslow \
			editfilter=1 \
			maxsites=100000 \
			overwrite=f \
			threads={threads} \
			outm={params.probedir_2mm}/"{{}}"_2mm.sam'

			find {params.probedir_nomm} -type f -name "*_nomm.sam" -empty | xargs --no-run-if-empty rm
			find {params.probedir_2mm} -type f -name "*_2mm.sam" -empty | xargs --no-run-if-empty rm

			for i in {params.probedir_nomm}/*_nomm.sam; do \
				[[ ! -f "$i" ]] && continue; \
				awk -F "\\t" '{{ print FILENAME, $0 }}' OFS="\\t" $i | \
				awk -F "/" '{{ print $NF }}' | \
				sed 's/_nomm.sam//1' | \
				cut -f1-2 >> {params.ispcrdir}/probes_hits_nomm_target.tsv; done

			for i in {params.probedir_2mm}/*_2mm.sam; do \
				[[ ! -f "$i" ]] && continue; \
				awk -F "\\t" '{{ print FILENAME, $0 }}' OFS="\\t" $i | \
				awk -F "/" '{{ print $NF }}' | \
				sed 's/_2mm.sam//1' | \
				cut -f1-2 >> {params.ispcrdir}/probes_hits_2mm_target.tsv; done

			rm -rf {params.ispcrdir}/probes	
		
		fi

		touch {output.status}
		"""

#Blast against per-genome amplicons, need to match primer name to probe name
#Blast against whole genome

rule probesearch_amplicon_offtarget:
	input:
		offtarget_genomes = rules.search_offtarget.output.status
	output:
		status = OUTDIR + "status/probesearch_offtarget.txt",
	conda: "../envs/ispcr.yaml"
	threads: 16
	params:
		ampdir = OUTDIR + "ispcr_offtarget/amplicon",
		probedir_nomm = directory(OUTDIR + "ispcr_offtarget/probes/no_mismatch"),
		probedir_2mm = directory(OUTDIR + "ispcr_offtarget/probes/2_mismatch"),
		ispcrdir = directory(OUTDIR + "ispcr_offtarget"),
		outdir = OUTDIR
	shell:
		"""
		if [[ -s {params.outdir}/genomes_offtarget_primerhits.txt ]]
		then
			cat {params.outdir}/genomes_offtarget_primerhits.txt | \
			parallel -j {threads} \
			'bbmap.sh \
			in={params.outdir}/probes.fasta \
			ref={params.ampdir}/"{{}}".fasta \
			nodisk \
			noheader=t \
			ambig=all \
			vslow \
			perfectmode \
			maxsites=100000 \
			overwrite=f \
			threads={threads} \
			outm={params.probedir_nomm}/"{{}}"_nomm.fasta'

			cat {params.outdir}/genomes_offtarget_primerhits.txt | \
			parallel -j {threads} \
			'bbmap.sh \
			in={params.outdir}/probes.fasta \
			ref={params.ampdir}/"{{}}".fasta \
			nodisk \
			noheader=t \
			ambig=all \
			vslow \
			editfilter=1 \
			maxsites=100000 \
			overwrite=f \
			threads={threads} \
			outm={params.probedir_2mm}/"{{}}"_2mm.fasta'

			find {params.probedir_nomm} -type f -name "*_nomm.fasta" -empty | xargs --no-run-if-empty rm
			find {params.probedir_2mm} -type f -name "*_2mm.fasta" -empty | xargs --no-run-if-empty rm

			for i in {params.probedir_nomm}/*_nomm.sam; do \
				[[ ! -f "$i" ]] && continue; \
				awk -F "\t" '{{ print FILENAME, $0 }}' OFS="\\t" $i | \
				awk -F "/" '{{ print $NF }}' | \
				sed 's/_nomm.sam//1' | \
				cut -f1-2 >> {params.ispcrdir}/probes_hits_nomm_offtarget.tsv; done
			

			for i in {params.probedir_2mm}/*_2mm.sam; do \
				[[ ! -f "$i" ]] && continue; \
				awk -F "\t" '{{ print FILENAME, $0 }}' OFS="\\t" $i | \
				awk -F "/" '{{ print $NF }}' | \
				sed 's/_2mm.sam//1' | \
				cut -f1-2 >> {params.ispcrdir}/probes_hits_2mm_offtarget.tsv; done
			
			rm -rf {params.ispcrdir}/probes

		fi

		touch {output.status}
		"""


rule collate_results:
	input:
		rules.probesearch_amplicon_target.output.status,
		rules.probesearch_amplicon_offtarget.output.status,
	output:
		status = OUTDIR + "status/collate.txt"
	params:
		outdir = OUTDIR,
		targetdir = OUTDIR + "ispcr_target",
		offtargetdir = OUTDIR + "ispcr_offtarget",
	shell:
		"""
		if [[ -s {params.outdir}/genomes_target_primerhits.txt ]]
		then
			cat {params.targetdir}/bed/*.bed > {params.targetdir}/target.bed
			cat {params.targetdir}/amplicon/*.fasta > {params.targetdir}/target_amplicons.fasta
			rm -rf {params.targetdir}/bed {params.targetdir}/amplicon
		fi

		if [[ -s {params.outdir}/genomes_offtarget_primerhits.txt ]]
		then
			cat {params.offtargetdir}/bed/*.bed > {params.offtargetdir}/offtarget.bed
			cat {params.offtargetdir}/amplicon/*.fasta > {params.offtargetdir}/offtarget_amplicons.fasta
			rm -rf {params.offtargetdir}/bed {params.offtargetdir}/amplicon
		fi

		touch {output.status}
		"""