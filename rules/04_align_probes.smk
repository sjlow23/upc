rule probe_blast_target:
	input:
		rules.collate_ispcr_target.output.status,
		probes = rules.prepare_probes.output.probes,
		amplicons = rules.collate_ispcr_target.output.targetdedup
	output:
		blast = OUTDIR + "ispcr_target/probes_blast.tsv",
		status = OUTDIR + "status/probesearch_target.txt",
	conda: "../envs/blast.yaml"
	threads: 16
	params:
	shell:
		"""
		sed -i "s/ /--/1" {input.amplicons}

		blastn -query {input.probes} \
			-subject {input.amplicons} \
			-task blastn-short \
			-qcov_hsp_perc 70 \
			-outfmt "6 qseqid sseqid slen length pident nident mismatch qseq sseq" > {output.blast}

		sed -i "s/--/\t/g" {output.blast}
		sed -i "s/:/\t/1" {output.blast}

		sed -i "s/--/ /1" {input.amplicons}
		
		touch {output.status}
		"""


rule probe_blast_offtarget:
	input:
		rules.collate_ispcr_offtarget.output.status,
		probes = OUTDIR + "probes.fasta",
		amplicons = rules.collate_ispcr_offtarget.output.offtargetdedup
	output:
		blast = OUTDIR + "ispcr_offtarget/probes_blast.tsv",
		status = OUTDIR + "status/probesearch_offtarget.txt",
	conda: "../envs/blast.yaml"
	threads: 16
	params:
	shell:
		"""
		sed -i "s/ /--/1" {input.amplicons}

		blastn -query {input.probes} \
			-subject {input.amplicons} \
			-task blastn-short \
			-qcov_hsp_perc 70 \
			-outfmt "6 qseqid sseqid slen length pident nident mismatch qseq sseq" > {output.blast}

		sed -i "s/--/\t/g" {output.blast}
		sed -i "s/:/\t/1" {output.blast}

		sed -i "s/--/ /1" {input.amplicons}

		touch {output.status}
		"""


checkpoint parse_blast:
	input:
		target = rules.probe_blast_target.output.blast,
		offtarget = rules.probe_blast_offtarget.output.blast,
	output:
		probe_target_dir = directory(OUTDIR + "ispcr_target/probes"),
		probe_offtarget_dir = directory(OUTDIR + "ispcr_offtarget/probes"),
		target = OUTDIR + "ispcr_target/probes_target.fasta",
		offtarget = OUTDIR + "ispcr_offtarget/probes_offtarget.fasta",
		status = OUTDIR + "status/parse_blast.txt"
	conda: "../envs/primer_mismatch.yaml"
	threads: 16
	params:
		probe_target = OUTDIR + "ispcr_target/probes",
		probe_offtarget = OUTDIR + "ispcr_offtarget/probes"
	shell:
		"""
		Rscript scripts/parse_blast.R {input.target} {input.offtarget} {output.target} {output.offtarget}

		seqkit split --by-id --id-regexp ".*--(.+)--.*" --by-id-prefix "" -O {params.probe_target} {output.target}
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
		mafft --thread {threads} {input.probes} > {output.aln}
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
		mafft --thread {threads} {input.probes} > {output.aln}
		"""
	
rule parse_probe_target:
	input:
		aln = OUTDIR + "ispcr_target/probe_alignments/{probe}.aln",
	output:
		tsv = OUTDIR + "ispcr_target/probe_parsed/{probe}.tsv"
	threads: 1
	conda: "../envs/primer_mismatch.yaml"
	shell:
		"""
		Rscript scripts/get_mismatches_probe.R {input.aln} {output.tsv}
		"""

rule parse_probe_offtarget:
	input:
		aln = OUTDIR + "ispcr_offtarget/probe_alignments/{probe}.aln",
	output:
		tsv = OUTDIR + "ispcr_offtarget/probe_parsed/{probe}.tsv"
	threads: 1
	conda: "../envs/primer_mismatch.yaml"
	shell:
		"""
		Rscript scripts/get_mismatches_probe.R {input.aln} {output.tsv}
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
		echo -e "probe_set\tori_probe\tprobe_seq\tgenome\tposition\treference\ttarget\tbinding_site\tref_sequence\tref_length\tmutation\ttype\n$(cat {output.target})" > {output.target}
		cat {output.target} | csvtk filter2 -t -f '$target!=0' > {output.target_clean}

		touch {output.status}
		"""


rule collate_probes_offtarget:
	input:
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
		mkdir -p {params.outdir}/offtarget

		cat {input.offtarget} | grep -v "reference" > {output.offtarget}
		echo -e "probe_set\tori_probe\tprobe_seq\tgenome\tposition\treference\ttarget\tbinding_site\tref_sequence\tref_length\tmutation\ttype\n$(cat {output.offtarget})" > {output.offtarget}
		cat {output.offtarget} | csvtk filter2 -t -f '$target!=0' > {output.offtarget_clean}

		touch {output.status}
		"""



rule summary_probes_target:
	input:
		target = rules.collate_probes_target.output.target,
	output:
		pergenome = OUTDIR + "stats_probes/target/stats_pergenome.tsv",
		percombo = OUTDIR + "stats_probes/target/stats_permutation_combo.tsv",
		perprobe = OUTDIR + "stats_probes/target/stats_perprobe.tsv",
		mismatchcount = OUTDIR + "stats_probes/target/stats_mismatches.tsv",
		grouped = OUTDIR + "stats_probes/target/stats_mismatches_grouped.tsv",
		status = OUTDIR + "status/summary_probes_target.txt"
	threads: 4
	conda: "../envs/primer_mismatch.yaml"
	params:
		targetdb = GENOMES_TARGET,
	shell:
		"""
		targetcount=$(ls {params.targetdb}/*.fna | wc -l)
		
		Rscript scripts/get_probe_summary.R {input.target} {output.pergenome} {output.percombo} {output.perprobe} {output.mismatchcount} {output.grouped} $targetcount
		touch {output.status}
		"""


rule summary_probes_offtarget:
	input:
		target = rules.collate_probes_offtarget.output.offtarget,
	output:
		pergenome = OUTDIR + "stats_probes/offtarget/stats_pergenome.tsv",
		percombo = OUTDIR + "stats_probes/offtarget/stats_permutation_combo.tsv",
		perprobe = OUTDIR + "stats_probes/offtarget/stats_perprobe.tsv",
		mismatchcount = OUTDIR + "stats_probes/offtarget/stats_mismatches.tsv",
		grouped = OUTDIR + "stats_probes/offtarget/stats_mismatches_grouped.tsv",
		status = OUTDIR + "status/summary_probes_offtarget.txt"
	threads: 4
	conda: "../envs/primer_mismatch.yaml"
	params:
		offtargetdb = GENOMES_OFFTARGET,
	shell:
		"""
		targetcount=$(ls {params.offtargetdb}/*.fna | wc -l)
		
		Rscript scripts/get_probe_summary.R {input.target} {output.pergenome} {output.percombo} {output.perprobe} {output.mismatchcount} {output.grouped} $targetcount
		touch {output.status}
		"""