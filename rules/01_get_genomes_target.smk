checkpoint download_target:
	output:
		genomedir = directory(GENOMES_TARGET),
	conda: "../envs/download.yaml"
	params:
		download_target = DOWNLOAD_TARGET,
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
	threads: 8
	shell:
		"""
		if [[ {params.download_target} == "yes" ]]
		then
			mkdir -p {params.metadir}
			
			if [[ {params.domain} == "viral" && {params.use_assembly} == "no" ]]
			then
				datasets download virus genome \
				taxon {params.spid} \
				--complete-only \
				--filename target.zip
				dataformat tsv virus-genome --force --package target.zip > {params.metadir}/metadata_target.tsv

				unzip target.zip -d {params.targetdir}
				rm target.zip

				# Subsample to max specified target sequences
				if [[ {params.max_target} != "no" ]]
				then
					cat {params.targetdir}/ncbi_dataset/data/genomic.fna | seqkit sample -n {params.max_target} -o {params.outdir}/genomic.fna
					seqkit split --by-id --by-id-prefix "" -O {params.targetdir} {params.outdir}/genomic.fna
					rm {params.outdir}/genomic.fna
				else
					seqkit split --by-id --by-id-prefix "" -O {params.targetdir} {params.targetdir}/ncbi_dataset/data/genomic.fna
				fi
				
				rm -rf {params.targetdir}/ncbi_dataset
				rm {params.targetdir}/README.md {params.targetdir}/md5sum.txt 
				cut -f1 {params.metadir}/metadata_target.tsv | grep -v Accession | awk '{{ print $0".fna" }}' | sort | uniq > {params.outdir}/target_genomes.txt
			else
				datasets download genome taxon {params.spid} \
				--assembly-source {params.db} \
				--assembly-version latest \
				--assembly-level {params.assembly_level} \
				--include genome,seq-report \
				--dehydrated \
				--filename target.zip
				
				unzip target.zip -d {params.targetdir}

				dataformat tsv genome --package target.zip > {params.metadir}/metadata_target_assembly.tsv
				datasets summary virus genome taxon {params.spid} --as-json-lines | \
					dataformat tsv virus-genome \
					--fields accession,geo-location,geo-region,isolate-collection-date,host-common-name,host-name,virus-common-name,virus-name,virus-tax-id \
						> {params.metadir}/metadata_target.tsv

				# Subsample if required
				if [[ {params.max_target} != "no" ]]
				then	
					shuf -n {params.max_target} {params.targetdir}/ncbi_dataset/fetch.txt > {params.targetdir}/ncbi_dataset/tmp
					mv {params.targetdir}/ncbi_dataset/tmp {params.targetdir}/ncbi_dataset/fetch.txt
				fi
				
				datasets rehydrate --directory {params.targetdir}
				rm target.zip

				for i in {params.targetdir}/ncbi_dataset/data/GC?_*/sequence_report.jsonl; do dataformat tsv genome-seq --inputfile $i | \
					awk -F "\\t" '{{ print $7, $1 }}' OFS="\\t" | \
					grep -v Accession >> {params.outdir}/assembly_accession.tmp; done
				sort -k1b,1 {params.outdir}/assembly_accession.tmp > {params.outdir}/assembly_accession.txt

				mv {params.targetdir}/ncbi_dataset/data/GC?_*/*.fna {params.targetdir}/
				rm -rf {params.targetdir}/ncbi_dataset {params.targetdir}/README.md {params.outdir}/assembly_accession.tmp

				find {params.targetdir} -type f -name 'GC*.fna' -exec bash -c 'mv "$1" "$(dirname "$1")/$(basename "$1" | cut -d"_" -f1,2).fna"' _ {{}} \;
				cut -f1 {params.metadir}/metadata_target_assembly.tsv | grep -v Accession | awk '{{ print $0".fna" }}' | sort | uniq > {params.outdir}/target_genomes.txt
			fi		
		else
			echo "Target genomes provided"
			cp -r {params.user_target} {params.targetdir}
		fi
		"""
