checkpoint split_amplicons_target:
	input:
		infile = rules.collate_ispcr_bbmap.output.status,
		target = rules.collate_ispcr_bbmap.output.targetdedup,
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
			--id-regexp '^\S+\s+(\w+)' --out-dir {output.alndir} {input.target}
		
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
		if [[ -s {input.offtarget} ]]
		then
			mkdir -p {params.alndir}
			seqkit split --by-id --by-id-prefix "" -j {threads} --two-pass --update-faidx \
				--id-regexp '^\S+\s+(\w+)' --out-dir {output.alndir} {input.offtarget}
		else
			mkdir -p {output.alndir}
		fi
		"""


rule align_split_target:
	input:
		fasta = OUTDIR + "ispcr_target/amplicons/{primer}.fasta",
	output:
		aln = OUTDIR + "ispcr_target/primer_alignments/{primer}.aln",
	threads: 4
	conda: "../envs/align.yaml"
	shell:
		"""
		ginsi --thread {threads} {input.fasta} > {output.aln}.tmp
		trimal -keepheader -in {output.aln}.tmp -out {output.aln} -gt 0.5
		rm {output.aln}.tmp
		"""


rule align_split_offtarget:
	input:
		fasta = OUTDIR + "ispcr_offtarget/amplicons/{primer}.fasta",
	output:
		aln = OUTDIR + "ispcr_offtarget/primer_alignments/{primer}.aln",
	threads: 4
	conda: "../envs/align.yaml"
	shell:
		"""
		ginsi --thread {threads} {input.fasta} > {output.aln}.tmp
		trimal -keepheader -in {output.aln}.tmp -out {output.aln} -gt 0.5
		rm {output.aln}.tmp
		"""

rule lookup_aln_target:
	input:
		target = rules.align_split_target.output.aln,
	output:
		tab = OUTDIR + "ispcr_target/primer_alignments/{primer}.tab",
	threads: 4
	conda: "../envs/seqkit.yaml"
	shell:
		"""
		if [[ -s {input.target} ]]; then
			seqkit fx2tab {input.target} | sed 's/ /\\t/g' | sed 's/\\t$//g' > {output.tab}
		else
			touch {output.tab}
		fi
		"""
		

## Account for no off-target hits
rule lookup_aln_offtarget:
	input:
		offtarget = rules.align_split_offtarget.output.aln,
	output:
		tab = OUTDIR + "ispcr_offtarget/primer_alignments/{primer}.tab",
	threads: 4
	conda: "../envs/seqkit.yaml"
	shell:
		"""
		if [[ -s {input.offtarget} ]]; then
			seqkit fx2tab {input.offtarget} | sed 's/ /\\t/g' | sed 's/\\t$//g' > {output.tab}
		else
			touch {output.tab}
		fi
		"""


rule merge_duplicates_target:
	input:
		aln = OUTDIR + "ispcr_target/primer_alignments/{primer}.tab",
	output:
		full = OUTDIR + "ispcr_target/primer_alignments/{primer}.tsv",
		aln = OUTDIR + "ispcr_target/primer_alignments/{primer}_all.aln",
	params:
		outdir = OUTDIR,
		lookup = OUTDIR + "ispcr_target/duplicates_target.tab",
		use_assembly = USE_ASSEMBLY
	conda: "../envs/primer_mismatch.yaml"
	shell:
		"""
		if [[ {params.use_assembly} == "yes" ]]
		then
			lookup={params.outdir}/target_assembly_accession.txt
		else
			lookup=NULL
		fi

		if [[ -s {params.lookup} ]]; then 
			Rscript scripts/merge_duplicates.R {params.lookup} {input.aln} {output.full} $lookup
			awk -F "," '{{ print ">"$1, $2, $3, $4, $5, "\\n"$6 }}' {output.full} > {output.aln}
		else
			Rscript scripts/merge_duplicates.R NULL {input.aln} {output.full} $lookup
			awk -F "," '{{ print ">"$1, $2, $3, $4, $5, "\\n"$6 }}' {output.full} > {output.aln}
		fi
		"""


rule merge_duplicates_offtarget:
	input:
		aln = OUTDIR + "ispcr_offtarget/primer_alignments/{primer}.tab",
	output:
		full = OUTDIR + "ispcr_offtarget/primer_alignments/{primer}.tsv",
		aln = OUTDIR + "ispcr_offtarget/primer_alignments/{primer}_all.aln",
	params:
		outdir = OUTDIR,
		lookup = OUTDIR + "ispcr_offtarget/duplicates_offtarget.tab",
		use_assembly = USE_ASSEMBLY
	conda: "../envs/primer_mismatch.yaml"
	shell:
		"""
		if [[ {params.use_assembly} == "yes" ]]
		then
			lookup={params.outdir}/offtarget_assembly_accession.txt
		else
			lookup=NULL
		fi

		if [[ -s {params.lookup} ]]; then 
			Rscript scripts/merge_duplicates.R {params.lookup} {input.aln} {output.full} $lookup
			awk -F "," '{{ print ">"$1, $2, $3, $4, $5, "\\n"$6 }}' {output.full} > {output.aln}
		else
			Rscript scripts/merge_duplicates.R NULL {input.aln} {output.full} $lookup
			awk -F "," '{{ print ">"$1, $2, $3, $4, $5, "\\n"$6 }}' {output.full} > {output.aln}
		fi
		"""


rule parse_aln_target:
	input:
		aln = OUTDIR + "ispcr_target/primer_alignments/{primer}_all.aln",
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
		aln = OUTDIR + "ispcr_offtarget/primer_alignments/{primer}_all.aln",
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
		target = OUTDIR + "stats_primers/target/summary_primers_target.tsv",
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


rule select_best_primermatch_target:
	input:
		target = rules.collate_primers_target.output.target,
	output:
		target = OUTDIR + "stats_primers/target/summary_primers_bestmatch.tsv",
		status = OUTDIR + "status/best_primer_target.txt"
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR + "stats_primers"
	shell:
		"""
		mkdir -p {params.outdir}/target
		Rscript scripts/get_best_match.R {input.target} {output.target} primer
		touch {output.status}
		"""


rule collate_primers_offtarget:
	input:
		offtarget = get_primers_offtarget_tsv,
		status = rules.collate_ispcr_offtarget.output.status,
	output:
		offtarget = OUTDIR + "stats_primers/offtarget/summary_primers_offtarget.tsv",
		offtarget_clean = OUTDIR + "summary_primers_offtarget_clean.tsv",
		status = OUTDIR + "status/collate_primers_offtarget.txt"
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR + "stats_primers"
	shell:
		"""
		if [[ -s {input.status} ]]
		then
			mkdir -p {params.outdir}/offtarget

			cat {input.offtarget} | grep -v "reference" > {output.offtarget}
			echo -e "primer_set\tori_primer\tprimer_seq\tgenome\tposition\treference\ttarget\tbinding_site\tref_sequence\tref_length\tmutation\ttype\n$(cat {output.offtarget})" > {output.offtarget}
			cat {output.offtarget} | csvtk filter2 -t -f '$target!=0' > {output.offtarget_clean}
		else
			touch {output.offtarget} {output.offtarget_clean}
		fi

		touch {output.status}
		"""


rule select_best_primermatch_offtarget:
	input:
		offtarget = rules.collate_primers_offtarget.output.offtarget,
	output:
		offtarget = OUTDIR + "stats_primers/offtarget/summary_primers_bestmatch.tsv",
		status = OUTDIR + "status/best_primer_offtarget.txt"
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR + "stats_primers"
	shell:
		"""
		mkdir -p {params.outdir}/offtarget

		if [[ -s {input.offtarget} ]]
		then
			Rscript scripts/get_best_match.R {input.offtarget} {output.offtarget} primer
		else
			touch {output.offtarget}
		fi

		touch {output.status}
		"""


rule summary_primers_target:
	input:
		target = rules.select_best_primermatch_target.output.target,
		#target = rules.collate_primers_target.output.target,
	output:
		pergenome = OUTDIR + "stats_primers/target/stats_pergenome.tsv",
		percombo = OUTDIR + "stats_primers/target/stats_permutation_combo.tsv",
		perprimer = OUTDIR + "stats_primers/target/stats_perprimer.tsv",
		mismatchcount = OUTDIR + "stats_primers/target/stats_mismatches.tsv",
		grouped = OUTDIR + "stats_primers/target/stats_mismatches_grouped.tsv",
		missing = OUTDIR + "stats_primers/target/stats_genomes_missing.tsv",
		genomestatus = OUTDIR + "stats_primers/target/genome_status.txt",
		status = OUTDIR + "status/summary_primers_target.txt"
	threads: 4
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR,
		targetdb = GENOMES_TARGET,
		subsample = SUBSAMPLE_TARGET,
	shell:
		"""
		targetcount=$(ls {params.targetdb}/*.fna | wc -l)
		
		if [[ {params.subsample} == "no" ]]
		then
			genomelist="target_genomes.txt"
		else
			genomelist="target_genomes_subsampled.txt"
		fi
	
		Rscript scripts/get_primer_summary.R {input.target} \
		{output.pergenome} \
		{output.percombo} \
		{output.perprimer} \
		{output.mismatchcount} \
		{output.grouped} \
		{output.missing} \
		{output.genomestatus} \
		{params.outdir}/$genomelist \
		$targetcount \
		"target"

		touch {output.status}
		"""


rule summary_primers_offtarget:
	input:
		status = rules.collate_ispcr_offtarget.output.status,
		#offtarget = rules.collate_primers_offtarget.output.offtarget,
		offtarget = rules.select_best_primermatch_offtarget.output.offtarget,
	output:
		pergenome = OUTDIR + "stats_primers/offtarget/stats_pergenome.tsv",
		percombo = OUTDIR + "stats_primers/offtarget/stats_permutation_combo.tsv",
		perprimer = OUTDIR + "stats_primers/offtarget/stats_perprimer.tsv",
		mismatchcount = OUTDIR + "stats_primers/offtarget/stats_mismatches.tsv",
		grouped = OUTDIR + "stats_primers/offtarget/stats_mismatches_grouped.tsv",
		missing = OUTDIR + "stats_primers/offtarget/stats_genomes_missing.tsv",
		genomestatus = OUTDIR + "stats_primers/offtarget/genome_status.txt",
		status = OUTDIR + "status/summary_primers_offtarget.txt"
	threads: 4
	conda: "../envs/primer_mismatch.yaml"
	params:
		offtargetdb = GENOMES_OFFTARGET,
		genomelist = OUTDIR + "offtarget_genomes.txt"
	shell:
		"""
		offtargetcount=$(wc -l {params.genomelist} | awk '{{ print $1 }}')
		
		if [[ -s {input.status} ]]
		then
			Rscript scripts/get_primer_summary.R {input.offtarget} \
			{output.pergenome} \
			{output.percombo} \
			{output.perprimer} \
			{output.mismatchcount} \
			{output.grouped} \
			{output.missing} \
			{output.genomestatus} \
			{params.genomelist} \
			$offtargetcount \
			"offtarget"
		else
			touch {output.pergenome} {output.percombo} {output.perprimer} {output.mismatchcount} {output.grouped} {output.missing} {output.genomestatus}
		fi

		touch {output.status}
		"""
