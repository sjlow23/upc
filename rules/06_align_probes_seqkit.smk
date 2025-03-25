rule locate_probes_target:
	input:
		status_collate = rules.collate_primers_target.output.status,
		probes = rules.prepare_probes.output.probes,
		amplicons = rules.collate_ispcr_bbmap.output.merged
	output:
		location = OUTDIR + "ispcr_target/probes_location.tsv",
		status = OUTDIR + "status/probesearch_target.txt",
	conda: "../envs/seqkit.yaml"
	threads: 16
	params:
		max_mismatch = PROBE_MAX_MISMATCH,
	shell:
		"""
		seqkit locate --max-mismatch {params.max_mismatch} --pattern-file {input.probes} {input.amplicons} > {output.location}	
		touch {output.status}
		"""

rule locate_probes_offtarget:
	input:
		status = rules.collate_ispcr_offtarget.output.status,
		status_collate = rules.collate_primers_offtarget.output.status,
		probes = rules.prepare_probes.output.probes,
		amplicons = rules.collate_ispcr_offtarget.output.offtargetamp
	output:
		location = OUTDIR + "ispcr_offtarget/probes_location.tsv",
		status = OUTDIR + "status/probesearch_offtarget.txt",
	conda: "../envs/seqkit.yaml"
	threads: 16
	params:
		max_mismatch = PROBE_MAX_MISMATCH,
	shell:
		"""
		#sed -i "s/ /--/1" {input.amplicons}
		seqkit locate --max-mismatch {params.max_mismatch} --pattern-file {input.probes} {input.amplicons} > {output.location}
		#sed -i "s/--/ /1" {input.amplicons}
		touch {output.status}
		"""

checkpoint parse_probes_target:
	input:
		target = rules.locate_probes_target.output.location
	output:
		probe_target_dir = directory(OUTDIR + "ispcr_target/probes"),
		target = OUTDIR + "ispcr_target/probes_target.fasta",
		status = OUTDIR + "status/parse_probes_target.txt"
	conda: "../envs/primer_mismatch.yaml"
	threads: 16
	params:
		probe_target = OUTDIR + "ispcr_target/probes",
	shell:
		"""
		Rscript scripts/parse_probes_locate.R {input.target} {output.target}
		seqkit split --by-id --id-regexp ".*--(.+)--.*" --by-id-prefix "" -O {params.probe_target} {output.target}
		touch {output.status}
		"""

checkpoint parse_probes_offtarget:
	input:
		offtarget = rules.locate_probes_offtarget.output.location
	output:
		probe_offtarget_dir = directory(OUTDIR + "ispcr_offtarget/probes"),
		offtarget = OUTDIR + "ispcr_target/probes_offtarget.fasta",
		status = OUTDIR + "status/parse_probes_offtarget.txt"
	conda: "../envs/primer_mismatch.yaml"
	threads: 16
	params:
		probe_offtarget = OUTDIR + "ispcr_offtarget/probes",
	shell:
		"""
		Rscript scripts/parse_probes_locate.R {input.offtarget} {output.offtarget}
		seqkit split --by-id --id-regexp ".*--(.+)--.*" --by-id-prefix "" -O {params.probe_offtarget} {output.offtarget}
		touch {output.status}
		"""

rule align_probes_target:
	input:
		probes = OUTDIR + "ispcr_target/probes/{probe}.fasta"
	output:
		aln = OUTDIR + "ispcr_target/probe_alignments/{probe}.aln"
	conda: "../envs/align.yaml"
	threads: 1
	params:
	shell:
		"""
		cp {input.probes} {output.aln}
		"""

rule align_probes_offtarget:
	input:
		probes = OUTDIR + "ispcr_offtarget/probes/{probe}.fasta"
	output:
		aln = OUTDIR + "ispcr_offtarget/probe_alignments/{probe}.aln"
	conda: "../envs/align.yaml"
	threads: 8
	params:
	shell:
		"""
		cp {input.probes} {output.aln}
		"""
	
rule parse_probe_target:
	input:
		aln = OUTDIR + "ispcr_target/probe_alignments/{probe}.aln",
	output:
		tsv = OUTDIR + "ispcr_target/probe_parsed/{probe}.tsv"
	threads: 1
	conda: "../envs/primer_mismatch.yaml"
	params:
		oriprobes = OUTDIR + "probes_expand.txt",
	shell:
		"""
		Rscript scripts/get_mismatches_probe.R {input.aln} {output.tsv} {params.oriprobes}
		"""

rule parse_probe_offtarget:
	input:
		aln = OUTDIR + "ispcr_offtarget/probe_alignments/{probe}.aln",
	output:
		tsv = OUTDIR + "ispcr_offtarget/probe_parsed/{probe}.tsv"
	threads: 1
	conda: "../envs/primer_mismatch.yaml"
	params:
		oriprobes = OUTDIR + "probes_expand.txt",
	shell:
		"""
		Rscript scripts/get_mismatches_probe.R {input.aln} {output.tsv} {params.oriprobes}
		"""

rule collate_probes_target:
	input:
		target = get_probes_tsv_t,
	output:
		target = OUTDIR + "stats_probes/target/summary_probes_target.tsv",
		target_clean = OUTDIR + "summary_probes_target_clean.tsv",
		status = OUTDIR + "status/collate_probes_target.txt"
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR + "stats_probes"
	shell:
		"""
		mkdir -p {params.outdir}/target

		cat {input.target} | grep -v "reference" > {output.target}
		echo -e "ori_probe\tori_primer\tprobe_seq\tgenome\tposition\treference\ttarget\tbinding_site\tref_sequence\tref_length\tmutation\ttype\n$(cat {output.target})" > {output.target}
		cat {output.target} | csvtk filter2 -t -f '$target!=0' > {output.target_clean}

		touch {output.status}
		"""

rule select_best_probe_target:
	input:
		target = rules.collate_probes_target.output.target,
	output:
		target = OUTDIR + "stats_probes/target/summary_probes_bestmatch.tsv",
		status = OUTDIR + "status/best_probe_target.txt"
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR + "stats_probes"
	shell:
		"""
		mkdir -p {params.outdir}/target
		Rscript scripts/get_best_match.R {input.target} {output.target} probe
		touch {output.status}
		"""

rule collate_probes_offtarget:
	input:
		offtarget_blast = rules.locate_probes_offtarget.output.location,
		offtarget = get_probes_tsv_ot,
	output:
		offtarget = OUTDIR + "stats_probes/offtarget/summary_probes_offtarget.tsv",
		offtarget_clean = OUTDIR + "summary_probes_offtarget_clean.tsv",
		status = OUTDIR + "status/collate_probes_offtarget.txt"
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR + "stats_probes"
	shell:
		"""
		if [[ -s {input.offtarget_blast} ]]
		then
			mkdir -p {params.outdir}/offtarget

			cat {input.offtarget} | grep -v "reference" > {output.offtarget}
			echo -e "ori_probe\tori_primer\tprobe_seq\tgenome\tposition\treference\ttarget\tbinding_site\tref_sequence\tref_length\tmutation\ttype\n$(cat {output.offtarget})" > {output.offtarget}
			cat {output.offtarget} | csvtk filter2 -t -f '$target!=0' > {output.offtarget_clean}
		else
			touch {output.offtarget} {output.offtarget_clean}
		fi

		touch {output.status}
		"""

rule select_best_probe_offtarget:
	input:
		offtarget = rules.collate_probes_offtarget.output.offtarget,
	output:
		offtarget = OUTDIR + "stats_probes/offtarget/summary_probes_bestmatch.tsv",
		status = OUTDIR + "status/best_probe_offtarget.txt"
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR + "stats_probes"
	shell:
		"""
		mkdir -p {params.outdir}/offtarget
		Rscript scripts/get_best_match.R {input.offtarget} {output.offtarget} probe
		touch {output.status}
		"""

rule summary_probes_target:
	input:
		target = rules.select_best_probe_target.output.target,
		#target = rules.collate_probes_target.output.target,
	output:
		pergenome = OUTDIR + "stats_probes/target/stats_pergenome.tsv",
		percombo = OUTDIR + "stats_probes/target/stats_permutation_combo.tsv",
		perprobe = OUTDIR + "stats_probes/target/stats_perprobe.tsv",
		mismatchcount = OUTDIR + "stats_probes/target/stats_mismatches.tsv",
		grouped = OUTDIR + "stats_probes/target/stats_mismatches_grouped.tsv",
		missing = OUTDIR + "stats_probes/target/stats_genomes_missing.tsv",
		genomestatus = OUTDIR + "stats_probes/target/genome_status.txt",
		status = OUTDIR + "status/summary_probes_target.txt"
	threads: 4
	conda: "../envs/primer_mismatch.yaml"
	params:
		outdir = OUTDIR,
		oriprobes = OUTDIR + "probes_expand.txt",
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

		Rscript scripts/get_probe_summary.R {input.target} \
		{output.pergenome} \
		{output.percombo} \
		{output.perprobe} \
		{output.mismatchcount} \
		{output.grouped} \
		{output.missing} \
		{output.genomestatus} \
		{params.outdir}/$genomelist \
		$targetcount \
		{params.oriprobes} \
		"target"

		touch {output.status}
		"""

rule summary_probes_offtarget:
	input:
		offtarget_blast = rules.locate_probes_offtarget.output.location,
		offtarget = rules.select_best_probe_offtarget.output.offtarget,
		#offtarget = rules.collate_probes_offtarget.output.offtarget,
	output:
		pergenome = OUTDIR + "stats_probes/offtarget/stats_pergenome.tsv",
		percombo = OUTDIR + "stats_probes/offtarget/stats_permutation_combo.tsv",
		perprobe = OUTDIR + "stats_probes/offtarget/stats_perprobe.tsv",
		mismatchcount = OUTDIR + "stats_probes/offtarget/stats_mismatches.tsv",
		grouped = OUTDIR + "stats_probes/offtarget/stats_mismatches_grouped.tsv",
		missing = OUTDIR + "stats_probes/offtarget/stats_genomes_missing.tsv",
		genomestatus = OUTDIR + "stats_probes/offtarget/genome_status.txt",
		status = OUTDIR + "status/summary_probes_offtarget.txt"
	threads: 4
	conda: "../envs/primer_mismatch.yaml"
	params:
		offtargetdb = GENOMES_OFFTARGET,
		genomelist = OUTDIR + "offtarget_genomes.txt",
		oriprobes = OUTDIR + "probes_expand.txt",
	shell:
		"""
		offtargetcount=$(wc -l {params.genomelist} | awk '{{ print $1 }}')
		
		if [[ -s {input.offtarget_blast} ]]
		then
			Rscript scripts/get_probe_summary.R {input.offtarget} \
			{output.pergenome} \
			{output.percombo} \
			{output.perprobe} \
			{output.mismatchcount} \
			{output.grouped} \
			{output.missing} \
			{output.genomestatus}
			{params.genomelist} \
			$offtargetcount \
			{params.oriprobes} \
			"offtarget"
		else
			touch {output.pergenome} {output.percombo} {output.perprobe} {output.mismatchcount} {output.grouped} {output.genomestatus}
		fi
		
		touch {output.status}
		"""