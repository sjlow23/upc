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
