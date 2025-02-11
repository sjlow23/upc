#https://stackoverflow.com/questions/67794112/snakemake-create-wildcards-from-output-directory-using-checkpoints

checkpoint split_amplicons_target:
	input:
		infile = rules.collate_ispcr_target.output.status,
		target = rules.collate_ispcr_target.output.targetdedup,
	output:
		alndir = directory(OUTDIR + "ispcr_target/amplicons"),
	conda: "../envs/seqkit.yaml"
	threads: 8
	params:
		outdir = OUTDIR,
		alndir = OUTDIR + "ispcr_target/amplicons"
	shell:
		"""
		mkdir -p {params.alndir}

		seqkit split --by-id --by-id-prefix "" -j {threads} --two-pass --update-faidx \
			--id-regexp '^\S+--(\w+)' --out-dir {output.alndir} {input.target}

		"""

checkpoint split_amplicons_offtarget:
	input:
		infile = rules.collate_ispcr_offtarget.output.status,
		offtarget = rules.collate_ispcr_offtarget.output.offtargetdedup,
	output:
		alndir = directory(OUTDIR + "ispcr_offtarget/amplicons"),
	conda: "../envs/seqkit.yaml"
	threads: 8
	params:
		outdir = OUTDIR,
		alndir = OUTDIR + "ispcr_offtarget/amplicons"	
	shell:
		"""
		mkdir -p {params.alndir}

		seqkit split --by-id --by-id-prefix "" -j {threads} --two-pass --update-faidx \
			--id-regexp '^\S+--(\w+)' --out-dir {output.alndir} {input.offtarget}
		"""

rule align_split_target:
	input:
		fasta = OUTDIR + "ispcr_target/amplicons/{primer}.fasta",
	output:
		aln = OUTDIR + "ispcr_target/primer_alignments/{primer}.aln",
	threads: 1
	conda: "../envs/align.yaml"
	shell:
		"""
		mafft --thread {threads} {input.fasta} > {output.aln}
		"""

rule align_split_offtarget:
	input:
		fasta = OUTDIR + "ispcr_offtarget/amplicons/{primer}.fasta",
	output:
		aln = OUTDIR + "ispcr_offtarget/primer_alignments/{primer}.aln",
	threads: 1
	conda: "../envs/align.yaml"
	shell:
		"""
		mafft --thread {threads} {input.fasta} > {output.aln}
		"""

rule parse_aln_target:
	input:
		aln = OUTDIR + "ispcr_target/primer_alignments/{primer}.aln",
	output:
		tsv = OUTDIR + "ispcr_target/primer_parsed/{primer}.tsv"
	threads: 1
	conda: "../envs/primer_mismatch.yaml"
	params:
		oriprimers = OUTDIR + "primers_expand.txt",
	shell:
		"""
		Rscript scripts/get_mismatches_primer.R {input.aln} {output.tsv} {params.oriprimers}
		"""

rule parse_aln_offtarget:
	input:
		aln = OUTDIR + "ispcr_offtarget/primer_alignments/{primer}.aln",
	output:
		tsv = OUTDIR + "ispcr_offtarget/primer_parsed/{primer}.tsv"
	threads: 1
	conda: "../envs/primer_mismatch.yaml"
	params:
		oriprimers = OUTDIR + "primers_expand.txt",
	shell:
		"""
		Rscript scripts/get_mismatches_primer.R {input.aln} {output.tsv} {params.oriprimers}
		"""



rule collate_primers_target:
	input:
		target = get_primers_target_tsv,
	output:
		target = OUTDIR + "summary_primers_target.tsv",
		target_clean = OUTDIR + "summary_primers_target_clean.tsv",
		status = OUTDIR + "status/collate_primers_target.txt"
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR + "stats_primers"
	shell:
		"""
		mkdir -p {params.outdir}/target

		cat {input.target} | grep -v "reference" > {output.target}
		echo -e "primer_set\tori_primer\tprimer_seq\tgenome\tposition\treference\ttarget\tbinding_site\tref_sequence\tref_length\tmutation\ttype\n$(cat {output.target})" > {output.target}
		cat {output.target} | csvtk filter2 -t -f '$target!=0' > {output.target_clean}

		touch {output.status}
		"""


rule collate_primers_offtarget:
	input:
		offtarget = get_primers_offtarget_tsv,
	output:
		offtarget = OUTDIR + "summary_primers_offtarget.tsv",
		offtarget_clean = OUTDIR + "summary_primers_offtarget_clean.tsv",
		status = OUTDIR + "status/collate_primers_offtarget.txt"
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR + "stats_primers"
	shell:
		"""
		mkdir -p {params.outdir}/offtarget

		cat {input.offtarget} | grep -v "reference" > {output.offtarget}
		echo -e "primer_set\tori_primer\tprimer_seq\tgenome\tposition\treference\ttarget\tbinding_site\tref_sequence\tref_length\tmutation\ttype\n$(cat {output.offtarget})" > {output.offtarget}
		cat {output.offtarget} | csvtk filter2 -t -f '$target!=0' > {output.offtarget_clean}

		touch {output.status}
		"""


rule summary_primers_target:
	input:
		target = rules.collate_primers_target.output.target,
	output:
		pergenome = OUTDIR + "stats_primers/target/stats_pergenome.tsv",
		percombo = OUTDIR + "stats_primers/target/stats_permutation_combo.tsv",
		perprimer = OUTDIR + "stats_primers/target/stats_perprimer.tsv",
		mismatchcount = OUTDIR + "stats_primers/target/stats_mismatches.tsv",
		grouped = OUTDIR + "stats_primers/target/stats_mismatches_grouped.tsv",
		status = OUTDIR + "status/summary_primers_target.txt"
	threads: 4
	conda: "../envs/primer_mismatch.yaml"
	params:
		targetdb = GENOMES_TARGET,
	shell:
		"""
		targetcount=$(ls {params.targetdb}/*.fna | wc -l)
		
		Rscript scripts/get_primer_summary.R {input.target} {output.pergenome} {output.percombo} {output.perprimer} {output.mismatchcount} {output.grouped} $targetcount
		touch {output.status}
		"""


rule summary_primers_offtarget:
	input:
		offtarget = rules.collate_primers_offtarget.output.offtarget,
	output:
		pergenome = OUTDIR + "stats_primers/offtarget/stats_pergenome.tsv",
		percombo = OUTDIR + "stats_primers/offtarget/stats_permutation_combo.tsv",
		perprimer = OUTDIR + "stats_primers/offtarget/stats_perprimer.tsv",
		mismatchcount = OUTDIR + "stats_primers/offtarget/stats_mismatches.tsv",
		grouped = OUTDIR + "stats_primers/offtarget/stats_mismatches_grouped.tsv",
		status = OUTDIR + "status/summary_primers_offtarget.txt"
	threads: 4
	conda: "../envs/primer_mismatch.yaml"
	params:
		offtargetdb = GENOMES_OFFTARGET,
	shell:
		"""
		offtargetcount=$(ls {params.offtargetdb}/*.fna | wc -l)
	
		Rscript scripts/get_primer_summary.R {input.offtarget} {output.pergenome} {output.percombo} {output.perprimer} {output.mismatchcount} {output.grouped} $offtargetcount
		touch {output.status}
		"""
