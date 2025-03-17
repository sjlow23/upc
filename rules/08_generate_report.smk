rule generate_report:
	input:
		primers_combo = rules.summary_primers_target.output.percombo,
        primers_missing = rules.summary_primers_target.output.missing,
		probes_combo = rules.summary_probes_target.output.percombo,
		probes_missing = rules.summary_probes_target.output.missing,
        no_amp = XXXXXXX,
	output:
		report = OUTDIR + "upc_report.html",
		status = OUTDIR + "status/generate_report.txt"
	conda: "../envs/report.yaml"
	threads: 16
	params:
	shell:
		"""
		Rscript scripts/print_html_report.R \
		{input.primers_combo} \
        {input.primers_missing} \
        {input.probes_combo} \
        {input.probes_missing} \
		{output.report} 
		
		touch {output.status}
		"""