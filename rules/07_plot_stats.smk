rule plot_stats:
	input:
		primers_combo = rules.summary_primers_target.output.percombo,
		primers_mutation = rules.summary_primers_target.output.perprimer,
		primers_status = rules.summary_primers_target.output.genomestatus,
		probes_combo = rules.summary_probes_target.output.percombo,
		probes_mutation = rules.summary_probes_target.output.perprobe,
		probes_grouped = rules.summary_probes_target.output.grouped,
		probes_status = rules.summary_probes_target.output.genomestatus,
	output:
		primerplot = OUTDIR + "stats_primers/primerplot.pdf",
		probeplot = OUTDIR + "stats_probes/probeplot.pdf",
		primer_combo_plot = OUTDIR + "primer_combo_plot.png",
		primer_mutation_plot = OUTDIR + "primer_mutation_plot.png",
		probe_combo_plot = OUTDIR + "probe_combo_plot.png",
		probe_mutation_plot = OUTDIR + "probe_mutation_plot.png",
		alluvial_plot = OUTDIR + "alluvial_plot.png",
		status = OUTDIR + "status/plot_stats.txt"
	conda: "../envs/plot.yaml"
	threads: 16
	params:
	shell:
		"""
		Rscript scripts/plot_stats.R \
		{input.primers_combo} \
		{input.primers_mutation} \
		{input.primers_status} \
		{input.probes_combo} \
		{input.probes_mutation} \
		{input.probes_grouped} \
		{input.probes_status} \
		{output.primerplot} \
		{output.probeplot} \
		{output.primer_combo_plot} \
		{output.primer_mutation_plot} \
		{output.probe_combo_plot} \
		{output.probe_mutation_plot} \
		{output.alluvial_plot}
		
		touch {output.status}
		"""

rule plot_piechart:
	input:
		primers = rules.summary_primers_target.output.genomestatus,
		probes = rules.summary_probes_target.output.genomestatus,
	output:
		primerpie = OUTDIR + "primer_piechart.png",
		probepie = OUTDIR + "probe_piechart.png",
		status = OUTDIR + "status/plot_piechart.txt"
	params:
		outdir = OUTDIR,
	conda: "../envs/plot.yaml"
	shell:
		"""
		Rscript scripts/plot_piechart.R {input.primers} {input.probes} {output.primerpie} {output.probepie}

		touch {output.status}
		"""

rule make_alignment_missing:
	input:
		rules.summary_primers_target.output.status,
		grouped = rules.summary_primers_target.output.grouped,
		primerstatus = rules.summary_primers_target.output.genomestatus,
		primers = OUTDIR + "primers_expand.txt"
	output:
		msaplot = OUTDIR + "msaplot_primers.png",
		status = OUTDIR + "status/make_alignment_missing.txt"
	params:
		outdir = OUTDIR,
		targetdir = GENOMES_TARGET,
		alndir = OUTDIR + "check_alignments",
		mode = "primer",
		segmentdir = OUTDIR + "split_segments",
		use_assembly = USE_ASSEMBLY
	conda: "../envs/msaplot.yaml",
	threads: CPU
	shell:
		"""
		if grep -qw "Not amplified" {input.primerstatus} && [[ {params.use_assembly} == "no" ]]; then
			mkdir -p {params.alndir}
			cut -f1 {input.primers} | sort | uniq > {params.outdir}/primerlist.txt

			if [[ {params.use_assembly} == "no" ]]
			then
				./scripts/missing_alignments.sh \
				{params.outdir}/primerlist.txt \
				{params.targetdir} \
				{params.alndir} \
				{input.grouped} \
				{input.primerstatus} \
				{input.primers} \
				{params.mode} \
				{threads}
			fi

				# mkdir -p {params.segmentdir}
				# ./scripts/split_fasta.sh {params.targetdir} {params.segmentdir} {params.outdir}/target_assembly_accession.txt
				# ./scripts/missing_alignments_assembly.sh \
				# {params.outdir}/primerlist.txt \
				# {params.targetdir} \
				# {params.alndir} \
				# {params.segmentdir} \
				# {input.grouped} \
				# {input.primerstatus} \
				# {input.primers} \
				# {params.mode} \
				# {threads}

			# Merge png plots
			magick convert {params.alndir}/*.png -bordercolor white -border 0x150 -append {output.msaplot}
	
			rm {params.outdir}/primerlist.txt {params.alndir}/*.png
		else
			touch {output.msaplot}
		fi

		touch {output.status}
		"""

rule make_probe_missing:
	input:
		rules.make_alignment_missing.output.status,
		rules.summary_probes_target.output.status,
		grouped = rules.summary_probes_target.output.grouped,
		probestatus = rules.summary_probes_target.output.genomestatus,
	output:
		msaplot = OUTDIR + "msaplot_probes.png",
		status = OUTDIR + "status/make_probe_missing.txt"
	params:
		outdir = OUTDIR,
		targetdir = GENOMES_TARGET,
		alndir = OUTDIR + "check_alignments",
		probes = OUTDIR + "probes_expand.txt",
		mode = "probe",
		use_assembly = USE_ASSEMBLY
	conda: "../envs/msaplot.yaml",
	threads: CPU
	shell:
		"""
		if grep -qw "Not present" {input.probestatus} && [[ {params.use_assembly} == "no" ]]; then
			mkdir -p {params.alndir}
			cut -f1 {params.probes} | sort | uniq > {params.outdir}/probelist.txt

			./scripts/missing_alignments.sh \
			{params.outdir}/probelist.txt \
			{params.targetdir} \
			{params.alndir} \
			{input.grouped} \
			{input.probestatus} \
			{params.probes} \
			{params.mode} \
			{threads}

			# Merge png plots
			magick convert {params.alndir}/*.png -bordercolor white -border 0x150 -append {output.msaplot}
		
			rm {params.outdir}/probelist.txt {params.alndir}/*.png
		else
			touch {output.msaplot}
		fi

		touch {output.status}
		"""