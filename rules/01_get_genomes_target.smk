checkpoint download_target:
	output:
		genomedir = directory(GENOMES_TARGET),
	conda: "../envs/download.yaml"
	params:
		download_target = DOWNLOAD_TARGET if "DOWNLOAD_TARGET" in globals() else [],
		user_target = USER_TARGET if "USER_TARGET" in globals() else [],
		spid = TARGET_SP_TAXID,
		targetdir = GENOMES_TARGET,
		domain = DOMAIN,
		db = DB,
		use_assembly = USE_ASSEMBLY,
		assembly_level = ASSEMBLY_LEVEL,
		outdir = OUTDIR,
		metadir = OUTDIR + "metadata/",
		max_target = SUBSAMPLE_TARGET,
		min_genome_size = MIN_GENOME_SIZE,
	threads: 8
	shell:
		"""
		if [[ {params.download_target} == "yes" ]]
		then
			mkdir -p {params.metadir} {params.outdir}
			
			./scripts/download_genomes_target.sh \
				{params.domain} \
				{params.use_assembly} \
				{params.spid} \
				{params.outdir} \
				{params.metadir} \
				{params.targetdir} \
				{params.max_target} \
				{params.min_genome_size} \
				{params.db} \
				{params.assembly_level}

		else
			mkdir -p {output.genomedir}
			cp {params.user_target}/* {output.genomedir}/
			
			for genome in {params.targetdir}/*.fna; do
				if [[ `seqkit stats $genome | tail -1 | awk '{{ print $7 }}' | sed 's/,//g'` -lt {params.min_genome_size} ]]
				then
					rm $genome
				fi
			done

			for i in {params.targetdir}/*.fna; do basename $i >> {params.outdir}/target_genomes.txt; done

			if [[ {params.max_target} != "no" ]]
			then
				shuf -n {params.max_target} {params.outdir}/target_genomes.txt > {params.outdir}/target_genomes_subsampled.txt
			fi
		fi
		"""

