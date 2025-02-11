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
		max_target = SUBSAMPLE_TARGET,
	threads: 8
	shell:
		"""
		if [[ {params.download_target} == "yes" ]]
		then
			if [[ {params.domain} == "viral" && {params.use_assembly} == "no" ]]
			then
				datasets download virus genome \
				taxon {params.spid} \
				--complete-only \
				--filename target.zip
				dataformat tsv virus-genome --force --package target.zip > {params.outdir}/metadata_target.tsv

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
			else
				datasets download genome taxon {params.spid} \
				--assembly-source {params.db} \
				--assembly-version latest \
				--assembly-level {params.assembly_level} \
				--include genome,seq-report \
				--dehydrated \
				--filename target.zip
				
				unzip target.zip -d {params.targetdir}

				dataformat tsv genome --package target.zip > {params.outdir}/metadata_target_assembly.tsv
				datasets summary virus genome taxon {params.spid} --as-json-lines | \
					dataformat tsv virus-genome \
					--fields accession,geo-location,geo-region,isolate-collection-date,host-common-name,host-name,virus-common-name,virus-name \
						> {params.outdir}/metadata_target.tsv

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
			fi
		
		elif [[ {params.download_target} == "no" ]]
		then
			echo "Target genomes provided"
			cp -r {params.user_target} {output.genomedir}
		else
			echo "No input found"
		fi
		
		"""


checkpoint download_offtarget:
	input:
	output:
		genomedir = directory(GENOMES_OFFTARGET),
	conda: "../envs/download.yaml"
	params:
		download_offtarget = DOWNLOAD_OFFTARGET,
		user_offtarget = USER_OFFTARGET if "USER_OFFTARGET" in globals() else [],
		spid = OFFTARGET_SP_TAXID,
		offtargetdir = GENOMES_OFFTARGET,
		domain = DOMAIN,
		db = DB,
		use_assembly = USE_ASSEMBLY,
		assembly_level = ASSEMBLY_LEVEL,
		outdir = OUTDIR,
		max_offtarget = SUBSAMPLE_OFFTARGET,
	threads: 4
	shell:
		"""
		if [[ {params.download_offtarget} == "yes" ]]
		then
			if [[ {params.domain} == "viral" && {params.use_assembly} == "no" ]]
			then
				datasets download virus genome taxon {params.spid} --complete-only --filename offtarget.zip
				dataformat tsv virus-genome --force --package offtarget.zip > {params.outdir}/metadata_offtarget.tsv

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

				dataformat tsv genome --package offtarget.zip > {params.outdir}/metadata_offtarget_assembly.tsv
				datasets summary virus genome taxon {params.spid} --as-json-lines | \
					dataformat tsv virus-genome --fields accession,geo-location,geo-region,isolate-collection-date,host-common-name,host-name,virus-name \
						> {params.outdir}/metadata_offtarget.tsv
				
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
			
		
		elif [[ {params.download_offtarget} == "no" ]]
		then
			echo "Off-target genomes provided"
			cp -r {params.user_offtarget} {output.genomedir}
		else
			echo "No input found"
		fi
		
		"""


rule prepare_primers:
	input:
		primers = PRIMERS,
	output:
		primers = OUTDIR + "primers.txt",
		primers_expand = OUTDIR + "primers_expand.txt",
		status = OUTDIR + "status/prepare_primers.txt",
	params:
		outdir = OUTDIR,
		probemode = PROBES,
		probes = OUTDIR + "probes.fasta"
	threads: 8
	conda: "../envs/download.yaml"
	shell:
		"""
		mkdir -p {params.outdir}
		
		awk -F "\\t" '{{ print $1, toupper($2), toupper($3) }}' OFS="\\t" {input.primers} > {params.outdir}/primerstmp.txt
		python scripts/expand_iupac.py {params.outdir}/primerstmp.txt {output.primers_expand} primers
		cut -f1,3-4 {output.primers_expand} > {output.primers}
		rm {params.outdir}/primerstmp.txt

		if [[ {params.probemode} == "yes" ]]
		then
			awk -F "\\t" '{{ print $1, toupper($4) }}' OFS="\\t" {input.primers} > {params.outdir}/probestmp.txt
			python scripts/expand_iupac.py {params.outdir}/probestmp.txt {params.outdir}/probes_expand.txt probes
			
			# Generate probes fasta
			cut -f1,3 probes_expand.txt | sort | uniq | awk 'BEGIN {{ prev=0 ; count=1 }} {{ if (prev==$1) count++; else {{ count=1;;prev=$1 }} print $1, $2, "probe"count }}' OFS="\\t" | \
				awk -F "\\t" '{{ print ">"$1"_"$3"\\n"$2 }}' > {params.probes}
			
			rm {params.outdir}/probestmp.txt
		fi
	
		touch {output.status}
		"""

rule prepare_probes:
	input:
		primers = PRIMERS,
	output:
		probes = OUTDIR + "probes.fasta",
		status = OUTDIR + "status/prepare_probes.txt",
	params:
		outdir = OUTDIR,
		probemode = PROBES,
	threads: 8
	conda: "../envs/download.yaml"
	shell:
		"""
		mkdir -p {params.outdir}
		
		awk -F "\\t" '{{ print $1, toupper($4) }}' OFS="\\t" {input.primers} > {params.outdir}/probestmp.txt
		python scripts/expand_iupac.py {params.outdir}/probestmp.txt {params.outdir}/probes_expand.txt probes
			
		# Generate probes fasta
		cut -f1,3 probes_expand.txt | sort | uniq | awk 'BEGIN {{ prev=0 ; count=1 }} {{ if (prev==$1) count++; else {{ count=1;;prev=$1 }} print $1, $2, "probe"count }}' OFS="\\t" | \
			awk -F "\\t" '{{ print ">"$1"_"$3"\\n"$2 }}' > {output.probes}
			
		rm {params.outdir}/probestmp.txt
	
		touch {output.status}
		"""