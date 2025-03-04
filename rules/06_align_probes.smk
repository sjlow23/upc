rule probe_blast_target:
	input:
		status_collate = rules.collate_primers_target.output.status,
		probes = rules.prepare_probes.output.probes,
		amplicons = rules.collate_ispcr_bbmap.output.merged
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
			-qcov_hsp_perc 50 \
			-outfmt "6 qseqid sseqid slen length pident nident mismatch sseq" > {output.blast}

		sed -i "s/--/\t/g" {output.blast}
		sed -i "s/:/\t/1" {output.blast}

		sed -i "s/--/ /1" {input.amplicons}
		
		touch {output.status}
		"""


rule probe_blast_offtarget:
	input:
		status = rules.collate_ispcr_offtarget.output.status,
		status_collate = rules.collate_primers_offtarget.output.status,
		probes = rules.prepare_probes.output.probes,
		amplicons = rules.collate_ispcr_offtarget.output.offtargetamp
	output:
		blast = OUTDIR + "ispcr_offtarget/probes_blast.tsv",
		status = OUTDIR + "status/probesearch_offtarget.txt",
	conda: "../envs/blast.yaml"
	threads: 16
	params:
	shell:
		"""
		if [[ -s {input.status} ]]
		then
			sed -i "s/ /--/1" {input.amplicons}
			blastn -query {input.probes} \
			-subject {input.amplicons} \
			-task blastn-short \
			-qcov_hsp_perc 50 \
			-outfmt "6 qseqid sseqid slen length pident nident mismatch sseq" > {output.blast}

			sed -i "s/--/\t/g" {output.blast}
			sed -i "s/:/\t/1" {output.blast}

			sed -i "s/--/ /1" {input.amplicons}
		else
			touch {output.blast}
		fi

		touch {output.status}
		
		"""


checkpoint parse_blast_target:
	input:
		target = rules.probe_blast_target.output.blast,
	output:
		probe_target_dir = directory(OUTDIR + "ispcr_target/probes"),
		target = OUTDIR + "ispcr_target/probes_target.fasta",
		status = OUTDIR + "status/parse_blast_target.txt"
	conda: "../envs/primer_mismatch.yaml"
	threads: 16
	params:
		probe_target = OUTDIR + "ispcr_target/probes",
		probes = OUTDIR + "probes.fasta"
	shell:
		"""
		Rscript scripts/parse_blast.R {params.probes} {input.target} {output.target}
		seqkit split --by-id --id-regexp ".*--(.+)--.*" --by-id-prefix "" -O {params.probe_target} {output.target}
		
		touch {output.status}
		"""


checkpoint parse_blast_offtarget:
	input:
		offtarget = rules.probe_blast_offtarget.output.blast,
	output:
		probe_offtarget_dir = directory(OUTDIR + "ispcr_offtarget/probes"),
		offtarget = OUTDIR + "ispcr_offtarget/probes_offtarget.fasta",
		status = OUTDIR + "status/parse_blast_offtarget.txt"
	conda: "../envs/primer_mismatch.yaml"
	threads: 16
	params:
		probe_offtarget = OUTDIR + "ispcr_offtarget/probes",
		probes = OUTDIR + "probes.fasta"
	shell:
		"""
		if [[ -s {input.offtarget} ]]
		then
			Rscript scripts/parse_blast.R {params.probes} {input.offtarget} {output.offtarget}
			seqkit split --by-id --id-regexp ".*--(.+)--.*" --by-id-prefix "" -O {output.probe_offtarget_dir} {output.offtarget}
		else
			mkdir -p {output.probe_offtarget_dir}
			touch {output.offtarget}
		fi

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
		ginsi --thread {threads} {input.probes} > {output.aln}
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
		ginsi --thread {threads} {input.probes} > {output.aln}
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
		echo -e "ori_probe\tprobe_seq\tgenome\tposition\treference\ttarget\tbinding_site\tref_sequence\tref_length\tmutation\ttype\n$(cat {output.target})" > {output.target}
		cat {output.target} | csvtk filter2 -t -f '$target!=0' > {output.target_clean}

		touch {output.status}
		"""


rule collate_probes_offtarget:
	input:
		offtarget_blast = rules.probe_blast_offtarget.output.blast,
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
			echo -e "ori_probe\tprobe_seq\tgenome\tposition\treference\ttarget\tbinding_site\tref_sequence\tref_length\tmutation\ttype\n$(cat {output.offtarget})" > {output.offtarget}
			cat {output.offtarget} | csvtk filter2 -t -f '$target!=0' > {output.offtarget_clean}
		else
			touch {output.offtarget} {output.offtarget_clean}
		fi

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
		missing = OUTDIR + "stats_probes/target/stats_genomes_missing.tsv",
		status = OUTDIR + "status/summary_probes_target.txt"
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

		Rscript scripts/get_probe_summary.R {input.target} \
		{output.pergenome} \
		{output.percombo} \
		{output.perprobe} \
		{output.mismatchcount} \
		{output.grouped} \
		{output.missing} \
		{params.outdir}/$genomelist \
		$targetcount \
		"target"

		touch {output.status}
		"""


rule summary_probes_offtarget:
	input:
		offtarget_blast = rules.probe_blast_offtarget.output.blast,
		offtarget = rules.collate_probes_offtarget.output.offtarget,
	output:
		pergenome = OUTDIR + "stats_probes/offtarget/stats_pergenome.tsv",
		percombo = OUTDIR + "stats_probes/offtarget/stats_permutation_combo.tsv",
		perprobe = OUTDIR + "stats_probes/offtarget/stats_perprobe.tsv",
		mismatchcount = OUTDIR + "stats_probes/offtarget/stats_mismatches.tsv",
		grouped = OUTDIR + "stats_probes/offtarget/stats_mismatches_grouped.tsv",
		missing = OUTDIR + "stats_probes/offtarget/stats_genomes_missing.tsv",
		status = OUTDIR + "status/summary_probes_offtarget.txt"
	threads: 4
	conda: "../envs/primer_mismatch.yaml"
	params:
		offtargetdb = GENOMES_OFFTARGET,
		genomelist = OUTDIR + "offtarget_genomes.txt"
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
			{params.genomelist} \
			$offtargetcount \
			"offtarget"
		else
			touch {output.pergenome} {output.percombo} {output.perprobe} {output.mismatchcount} {output.grouped}
		fi
		
		touch {output.status}
		"""