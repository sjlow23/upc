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
			
			if [[ {params.domain} == "viral" && {params.use_assembly} == "no" ]]
			then
				datasets download virus genome \
				taxon {params.spid} \
				--complete-only \
				--filename {params.outdir}/target.zip
				dataformat tsv virus-genome --force --package {params.outdir}/target.zip > {params.metadir}/metadata_target.tsv

				unzip {params.outdir}/target.zip -d {params.targetdir}
				rm {params.outdir}/target.zip

				# Subsample to max specified target sequences
				if [[ {params.max_target} != "no" ]]
				then
					cat {params.targetdir}/ncbi_dataset/data/genomic.fna | \
						seqkit seq --min-len {params.min_genome_size} | \
						seqkit sample -n {params.max_target} -o {params.outdir}/genomic.fna
					seqkit split --by-id --by-id-prefix "" -O {params.targetdir} {params.outdir}/genomic.fna
					rm {params.outdir}/genomic.fna
					ls {params.targetdir}/*.fna | awk -F "/" '{{ print $NF }}' > {params.outdir}/target_genomes_subsampled.txt
				else
					seqkit split --by-id --by-id-prefix "" -O {params.targetdir} {params.targetdir}/ncbi_dataset/data/genomic.fna
				fi
				
				rm -rf {params.targetdir}/ncbi_dataset
				rm {params.targetdir}/README.md {params.targetdir}/md5sum.txt 

				for genome in {params.targetdir}/*.fna; do
					if [[ `seqkit stats $genome | tail -1 | awk '{{ print $7 }}' | sed 's/,//g'` -lt {params.min_genome_size} ]]
					then
						rm $genome
					fi
				done

				awk -F "\\t" '$27 >= {params.min_genome_size}' {params.metadir}/metadata_target.tsv | \
					cut -f1 | \
					grep -v Accession | \
					awk '{{ print $0".fna" }}' | \
					sort | uniq > {params.outdir}/target_genomes.txt
			else
				datasets download genome taxon {params.spid} \
				--assembly-source {params.db} \
				--assembly-version latest \
				--assembly-level {params.assembly_level} \
				--include genome,seq-report \
				--dehydrated \
				--filename {params.outdir}/target.zip
				
				unzip {params.outdir}/target.zip -d {params.targetdir}

				dataformat tsv genome --package {params.outdir}/target.zip > {params.metadir}/metadata_target_assembly.tsv
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
				rm {params.outdir}/target.zip

				for i in {params.targetdir}/ncbi_dataset/data/GC?_*/sequence_report.jsonl; do dataformat tsv genome-seq --inputfile $i | \
					awk -F "\\t" '{{ print $7, $1 }}' OFS="\\t" | \
					grep -v Accession >> {params.outdir}/assembly_accession.tmp; done
				sort -k1b,1 {params.outdir}/assembly_accession.tmp > {params.outdir}/assembly_accession.txt

				mv {params.targetdir}/ncbi_dataset/data/GC?_*/*.fna {params.targetdir}/
				rm -rf {params.targetdir}/ncbi_dataset {params.targetdir}/README.md {params.outdir}/assembly_accession.tmp

				find {params.targetdir} -type f -name 'GC*.fna' -exec bash -c 'mv "$1" "$(dirname "$1")/$(basename "$1" | cut -d"_" -f1,2).fna"' _ {{}} \;

				for genome in {params.targetdir}/*.fna; do
					if [[ `seqkit stats $genome | tail -1 | awk '{{ print $7 }}' | sed 's/,//g'` -lt {params.min_genome_size} ]]
					then
						rm $genome
					fi
				done

				awk -F "\\t" '$27 >= {params.min_genome_size}' {params.metadir}/metadata_target_assembly.tsv | \
					cut -f1 | \
					grep -v Accession | \
					awk '{{ print $0".fna" }}' | \
					sort | uniq > {params.outdir}/target_genomes.txt
			fi		
		else
			echo "Target genomes provided"
			mkdir -p {output.genomedir}
			cp {params.user_target}/* {output.genomedir}/
			
			for genome in {params.targetdir}/*.fna; do
				if [[ `seqkit stats $genome | tail -1 | awk '{{ print $7 }}' | sed 's/,//g'` -lt {params.min_genome_size} ]]
				then
					rm $genome
				fi
			done

			for i in {params.targetdir}/*.fna; do basename $i >> {params.outdir}/target_genomes.txt; done
		fi
		"""

