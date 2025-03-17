rule plot_stats:
	input:
		primers_combo = rules.summary_primers_target.output.percombo,
		primers_mutation = rules.summary_primers_target.output.perprimer,
		probes_combo = rules.summary_probes_target.output.percombo,
		probes_mutation = rules.summary_probes_target.output.perprobe,
	output:
		primerplot = OUTDIR + "stats_primers/primerplot.pdf",
		probeplot = OUTDIR + "stats_probes/probeplot.pdf",
		status = OUTDIR + "status/plot_stats.txt"
	conda: "../envs/plot.yaml"
	threads: 16
	params:
	shell:
		"""
		Rscript scripts/plot_stats.R \
		{input.primers_combo} \
		{input.primers_mutation} \
		{input.probes_combo} \
		{input.probes_mutation} \
		{output.primerplot} \
		{output.probeplot}
		
		touch {output.status}
		"""


rule summarise_genomes_primers:
	input:
		rules.plot_stats.output.status,
		mismatch = OUTDIR + "stats_primers/target/stats_mismatches.tsv",
		missing = OUTDIR + "stats_primers/target/stats_genomes_missing.tsv",
		primers = OUTDIR + "primers_expand.txt"
	output:
		stats = OUTDIR + "stats_primers/target/genome_status.txt",
		status = OUTDIR + "status/summarise_genomes_primers.txt"
	params:
		outdir = OUTDIR,
		statsdir = OUTDIR + "stats_primers/target",
		subsample = SUBSAMPLE_TARGET,
	conda: "../envs/primer_mismatch.yaml"
	shell:
		"""
		if [[ {params.subsample} != "no" ]]; then
			genomelist="{params.outdir}/target_genomes_subsampled.txt"
		else
			genomelist="{params.outdir}/target_genomes.txt"
		fi

		Rscript scripts/get_genome_summary.R {input.mismatch} {input.missing} {input.primers} {output.stats} $genomelist primers

		touch {output.status}
		"""


rule summarise_genomes_probes:
	input:
		rules.plot_stats.output.status,
		mismatch = OUTDIR + "stats_probes/target/stats_mismatches.tsv",
		missing = OUTDIR + "stats_probes/target/stats_genomes_missing.tsv",
	output:
		stats = OUTDIR + "stats_probes/target/genome_status.txt",
		status = OUTDIR + "status/summarise_genomes_probes.txt"
	params:
		outdir = OUTDIR,
		statsdir = OUTDIR + "stats_probes/target",
		subsample = SUBSAMPLE_TARGET,
		probes = OUTDIR + "probes_expand.txt"
	conda: "../envs/primer_mismatch.yaml"
	shell:
		"""
		if [[ {params.subsample} != "no" ]]; then
			genomelist="{params.outdir}/target_genomes_subsampled.txt"
		else
			genomelist="{params.outdir}/target_genomes.txt"
		fi

		Rscript scripts/get_genome_summary.R {input.mismatch} {input.missing} {params.probes} {output.stats} $genomelist probes

		touch {output.status}
		"""



rule make_alignment_missing:
	input:
		rules.summarise_genomes_primers.output.status,
		primers = rules.summarise_genomes_primers.output.stats,
		probes = rules.summarise_genomes_probes.output.stats,
		primerlist = OUTDIR + "primers_expand.txt"
	output:
		status = OUTDIR + "status/make_alignment_missing.txt"
	params:
		outdir = OUTDIR,
		targetdir = GENOMES_TARGET,
		alndir = OUTDIR + "check_alignments"
	conda: "../envs/align.yaml",
	threads: CPU
	shell:
		"""
		mkdir -p {params.alndir}

		cut -f1 {input.primerlist} | sort | uniq > {params.outdir}/primerlist.txt

		while read primer
		do
			grep "not amplified" {input.primers} | grep -w "$primer" | cut -f1 > {params.outdir}/keep_"$primer".txt
			grep "not present" {input.probes} | grep -w "$primer" | cut -f1 >> {params.outdir}/keep_"$primer".txt

			if [[ -s {params.outdir}/keep_"$primer".txt ]]; then
				grep -i perfect {input.primers} | shuf -n20 | cut -f1 >> {params.outdir}/keep_"$primer".txt
				grep -i mismatch {input.primers} | shuf -n10 | cut -f1 >> {params.outdir}/keep_"$primer".txt

				sort {params.outdir}/keep_"$primer".txt | uniq > {params.outdir}/keep2_"$primer".txt

				for i in `cat {params.outdir}/keep2_"$primer".txt`; do \
					cat {params.targetdir}/"$i".fna >> {params.alndir}/check_"$primer".fasta; done
				einsi --thread {threads} {params.alndir}/check_"$primer".fasta > {params.alndir}/check_"$primer".aln
				rm {params.outdir}/keep_"$primer".txt {params.outdir}/keep2_"$primer".txt
			fi
		done < {params.outdir}/primerlist.txt

		touch {output.status}
		"""
