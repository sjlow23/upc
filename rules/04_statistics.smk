rule primer_coverage_target:
	input:
		rules.collate_results.output.status,
		tree = OUTDIR + "phylogeny/combined_target.nwk",
		primers = rules.prepare_primers.output.primers_expand,
		phylo = OUTDIR + "status/phylogeny.txt"
	output:
		primerinfo = OUTDIR + "primer_statistics/primer_info.tsv",
		fullstats = OUTDIR + "primer_statistics/all_statistics.tsv",
		primerstats = OUTDIR + "primer_statistics/primer_statistics.tsv",
		genomestats = OUTDIR + "primer_statistics/genome_statistics.tsv",
		barplot = OUTDIR + "primer_statistics/barplot.pdf",
		treeplot = OUTDIR + "primer_statistics/treeplot.pdf",
		status = OUTDIR + "status/primerstats.txt"
	conda: "../envs/stats.yaml"
	threads: 1
	params:
		outdir = OUTDIR,
		bed = OUTDIR + "ispcr_target/target.bed",
		fasta = OUTDIR + "ispcr_target/target_amplicons.fasta",
		ispcr_dir = OUTDIR + "ispcr_target",
		output_dir = OUTDIR + "primer_statistics",
		metadata = OUTDIR + "phylogeny/metadata_subset.tsv"
	shell:
		"""
		grep ">" {params.fasta} | sed "s/>//g" | sed "s/ /\\t/g" | sed "s/:/\\t/1" > {output.primerinfo}
		awk -F "\\t" '{{$2=gensub(/[+-]/,"\\t","g",$2)}}1' OFS="\\t" {output.primerinfo} > {output.primerinfo}.2
		mv {output.primerinfo}.2 {output.primerinfo}

		if [[ -s {params.metadata} ]]
		then
			metadata={params.metadata}
		else
			metadata=NULL
		fi

		if [[ -s {params.outdir}/assembly_accession.txt ]]
		then
			accession={params.outdir}/assembly_accession.txt
		else
			accession=NULL
		fi

		Rscript scripts/parse_new.R {output.primerinfo} {params.outdir}/genomes_target.txt \
			{params.bed} \
			{input.tree} \
			$metadata \
			$accession \
			{input.primers} \
			{output.fullstats} \
			{output.primerstats} \
			{output.genomestats} \
			{output.barplot} \
			{output.treeplot} 
	
		touch {output.status}
		"""