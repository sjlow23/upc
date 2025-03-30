checkpoint download_offtarget:
	input:
	output:
		genomedir = directory(GENOMES_OFFTARGET_TMP),
	conda: "../envs/download.yaml"
	params:
		download_offtarget = DOWNLOAD_OFFTARGET if "DOWNLOAD_OFFTARGET" in globals() else [],
		user_offtarget = USER_OFFTARGET if "USER_OFFTARGET" in globals() else [],
		spid = OFFTARGET_SP_TAXID,
		offtargetdir = directory(GENOMES_OFFTARGET_TMP),
		domain = DOMAIN,
		db = DB,
		use_assembly = USE_ASSEMBLY,
		assembly_level = ASSEMBLY_LEVEL,
		outdir = OUTDIR,
		metadir = OUTDIR + "metadata/",
		max_offtarget = SUBSAMPLE_OFFTARGET,
	threads: 4
	shell:
		"""
		if [[ {params.download_offtarget} == "yes" ]]
		then
			mkdir -p {params.metadir} {params.outdir}
			
			./scripts/download_genomes_offtarget.sh \
				{params.domain} \
				{params.use_assembly} \
				{params.spid} \
				{params.outdir} \
				{params.metadir} \
				{params.offtargetdir} \
				{params.max_offtarget} \
				{params.db} \
				{params.assembly_level}

		else 
			mkdir -p {output.genomedir}
			cp {params.user_offtarget}/* {output.genomedir}/
			for i in {params.offtargetdir}/*.fna; do basename $i >> {params.outdir}/offtarget_genomes.txt; done

			if [[ {params.max_offtarget} != "no" ]]
			then
				shuf -n {params.max_offtarget} {params.outdir}/offtarget_genomes.txt > {params.outdir}/offtarget_genomes_subsampled.txt
			fi
		fi
		"""


