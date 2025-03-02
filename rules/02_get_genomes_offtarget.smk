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
			mkdir -p {params.metadir}
			if [[ {params.domain} == "viral" && {params.use_assembly} == "no" ]]
			then
				datasets download virus genome taxon {params.spid} --complete-only --filename offtarget.zip
				dataformat tsv virus-genome --force --package offtarget.zip > {params.metadir}/metadata_offtarget.tsv

				unzip offtarget.zip -d {params.offtargetdir}
				rm offtarget.zip

				# Subsample to max specified offtarget sequences
				if [[ {params.max_offtarget} != "no" ]]
				then
					cat {params.offtargetdir}/ncbi_dataset/data/genomic.fna | seqkit sample -n {params.max_offtarget} -o {params.outdir}/genomic.fna
					seqkit split --by-id --by-id-prefix "" -O {params.offtargetdir} {params.outdir}/genomic.fna
					rm {params.outdir}/genomic.fna
				else
					seqkit split --by-id --by-id-prefix "" -O {params.offtargetdir} {params.offtargetdir}/ncbi_dataset/data/genomic.fna
				fi

				rm -rf {params.offtargetdir}/ncbi_dataset
				rm {params.offtargetdir}/README.md {params.offtargetdir}/md5sum.txt 

			else
				datasets download genome taxon {params.spid} \
				--assembly-source {params.db} \
				--assembly-version latest \
				--assembly-level {params.assembly_level} \
				--include genome \
				--dehydrated \
				--filename offtarget.zip
				
				unzip offtarget.zip -d {params.offtargetdir}

				dataformat tsv genome --package offtarget.zip > {params.metadir}/metadata_offtarget_assembly.tsv
				datasets summary virus genome taxon {params.spid} --as-json-lines | \
					dataformat tsv virus-genome --fields accession,geo-location,geo-region,isolate-collection-date,host-common-name,host-name,virus-name \
						> {params.metadir}/metadata_offtarget.tsv
				
				# Subsample if required
				if [[ {params.max_offtarget} != "no" ]]
				then
					shuf -n {params.max_offtarget} {params.offtargetdir}/ncbi_dataset/fetch.txt > {params.offtargetdir}/ncbi_dataset/tmp
					mv {params.offtargetdir}/ncbi_dataset/tmp {params.offtargetdir}/ncbi_dataset/fetch.txt
				fi

				datasets rehydrate --directory {params.offtargetdir}
				
				rm offtarget.zip
				mv {params.offtargetdir}/ncbi_dataset/data/GC?_*/*.fna {params.offtargetdir}/
				rm -rf {params.offtargetdir}/ncbi_dataset {params.offtargetdir}/README.md

				find {params.offtargetdir} -type f -name 'GC*.fna' -exec bash -c 'mv "$1" "$(dirname "$1")/$(basename "$1" | cut -d"_" -f1,2).fna"' _ {{}} \;

			fi
			find {params.offtargetdir} -name "*.fna" | awk -F "/" '{{ print $NF }}' >> {params.outdir}/offtarget_genomes.txt

		else 
			echo "Off-target genomes provided"
			mkdir -p {output.genomedir}
			cp {params.user_offtarget}/* {output.genomedir}/
			for i in {params.offtargetdir}/*.fna; do basename $i >> {params.outdir}/offtarget_genomes.txt; done
		fi
		"""


