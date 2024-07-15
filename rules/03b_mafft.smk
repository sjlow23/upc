rule infer_phylogeny:
	input:
		inputdir = GENOMES_TARGET,
		genomelist = OUTDIR + "genomes_target.txt"
	output:
		aln = OUTDIR + "phylogeny/combined_target.aln",
		tree = OUTDIR + "phylogeny/combined_target.nwk",
		status = OUTDIR + "status/phylogeny.txt"
	conda: "../envs/align.yaml"
	params:
		outdir = OUTDIR + "phylogeny",
		derep = DEREP,
	threads: 16
	shell:
		"""
		# Set randomly selected target sequence as reference
		ref=$(find {input.inputdir}/*.fna | shuf -n1)
		
		if [[ {params.derep} == "no" ]]
		then
			cat {input.inputdir}/*.fna > {params.outdir}/combined_target.fasta
			mafft --thread {threads} --6merpair --keeplength --addfragments {params.outdir}/combined_target.fasta $ref > {output.aln}
		elif [[ {params.derep} == "yes" ]]
		then
			seqkit rmdup -s {params.outdir}/combined_target_nodups.aln > {output.aln}
			rm {params.outdir}/combined_target_nodups.aln
		fi

		export OMP_NUM_THREADS={threads}
		FastTreeMP -gamma -nt -gtr {output.aln} > {output.tree}
		touch {output.status}

		"""


rule tree_metadata:
	input:
		infile = rules.infer_phylogeny.output.status,
	output:
		subset_meta = OUTDIR + "phylogeny/metadata_subset.tsv"
	params:
		metadata = OUTDIR + "metadata_target.tsv"
	shell:
		"""
		cp {params.metadata} {output.subset_meta}
		
		"""
