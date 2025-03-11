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


rule summarise_genomes:
	input:
		rules.plot_stats.output.status,
		mismatch = OUTDIR + "stats_primers/target/stats_mismatches.tsv",
		missing = OUTDIR + "stats_primers/target/stats_genomes_missing.tsv",
	output:
		stats = OUTDIR + "stats_primers/target/genome_status.txt",
		status = OUTDIR + "status/summarise_genomes.txt"
	params:
		outdir = OUTDIR,
		statsdir = OUTDIR + "stats_primers/target",
		subsample = SUBSAMPLE_TARGET,
	shell:
		"""
		grep -v genome {input.mismatch} | cut -f1 | sort | uniq | awk '{{ print $0, "with mismatch" }}' OFS="\\t" > {output.stats}
		grep -v missing {input.missing} | cut -f5 | sed 's/, /\\n/g' | grep -v "^$" | sort | uniq | awk '{{ print $0, "not amplified" }}' OFS="\\t" >> {output.stats}

		cut -f1 {output.stats} > {params.statsdir}/tmp
		
		if [[ {params.subsample} == "yes" ]]; then
			sed 's/.fna//g' {params.outdir}/target_genomes_subsampled.txt > {params.statsdir}/genomestmp
		else
			sed 's/.fna//g' {params.outdir}/target_genomes.txt > {params.statsdir}/genomestmp; fi


		if [[ `wc -l {params.statsdir}/tmp | awk '{{ print $1 }}'` != `wc -l {params.statsdir}/genomestmp | awk '{{ print $1 }}'` ]]
		then
			grep -vwF -f {params.statsdir}/tmp {params.statsdir}/genomestmp | awk '{{ print $0, "perfect match" }}' OFS="\\t" >> {output.stats}
		else
			continue 
		fi
				
		rm {params.statsdir}/tmp {params.statsdir}/genomestmp
		touch {output.status}

		"""