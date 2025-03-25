rule generate_report:
	input:
		primerlist = OUTDIR + "primers_expand.txt",
		oriprimers = PRIMERS,
		primers_combo = rules.summary_primers_target.output.percombo,
        primers_missing = rules.summary_primers_target.output.missing,
		probes_combo = rules.summary_probes_target.output.percombo,
		probes_missing = rules.summary_probes_target.output.missing,
        primer_combo_plot  = rules.plot_stats.output.primer_combo_plot,
		primer_mutation_plot = rules.plot_stats.output.primer_mutation_plot,
		probe_combo_plot = rules.plot_stats.output.probe_combo_plot,
		probe_mutation_plot = rules.plot_stats.output.probe_mutation_plot,
		primerpie = rules.plot_piechart.output.primerpie,
		probepie = rules.plot_piechart.output.probepie,
		alluvial = rules.plot_stats.output.alluvial_plot,
		primermsaplot = rules.make_alignment_missing.output.msaplot,
		probemsaplot = rules.make_probe_missing.output.msaplot,
	output:
		report = OUTDIR + "upc_report.html",
		status = OUTDIR + "status/generate_report.txt"
	conda: "../envs/report.yaml"
	params:
		probelist = OUTDIR + "probes_expand.txt"
	threads: 16
	shell:
		"""
		Rscript scripts/print_html_report.R \
		{input.primerlist}	\
		{params.probelist} \
		{input.primers_combo} \
        {input.primers_missing} \
        {input.probes_combo} \
        {input.probes_missing} \
		{input.primer_combo_plot} \
		{input.primer_mutation_plot} \
		{input.probe_combo_plot} \
		{input.probe_mutation_plot} \
		{input.primerpie} \
		{input.probepie} \
		{input.alluvial} \
		{input.primermsaplot} \
		{input.probemsaplot} \
		{input.oriprimers} \
		{output.report} 
		
		rm {input.primer_combo_plot} {input.primer_mutation_plot}
		rm {input.probe_combo_plot} {input.probe_mutation_plot}
		rm {input.primerpie} {input.probepie}
		rm {input.alluvial} {input.primermsaplot} {input.probemsaplot}

		touch {output.status}
		"""