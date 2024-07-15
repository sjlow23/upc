checkpoint download_target:
	input:
		infile = "myseq.fasta",
	output:
		genomedir = directory(GENOMES_TARGET),
		genomelist = OUTDIR + "genomes_target.txt",
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
	threads: 4
	shell:
		"""
		if [[ {params.download_target} == "yes" ]]
		then
			if [[ {params.domain} == "viral" && {params.use_assembly} == "no" ]]
			then
				datasets download virus genome taxon {params.spid} --complete-only --filename target.zip
				dataformat tsv virus-genome --package target.zip > {params.outdir}/metadata_target.tsv

				unzip target.zip -d {params.targetdir}
				rm target.zip
				seqkit split --by-id -O {params.targetdir} {params.targetdir}/ncbi_dataset/data/genomic.fna
				ls {params.targetdir}/*.fna | parallel -j {threads} 'rename "genomic.part_" "" "{{}}"'  
				rm -rf {params.targetdir}/ncbi_dataset
				rm {params.targetdir}/README.md
			else
				datasets download genome taxon {params.spid} \
				--assembly-source {params.db} \
				--assembly-version latest \
				--assembly-level {params.assembly_level} \
				--include genome \
				--dehydrated \
				--filename target.zip
				dataformat tsv genome --package target.zip > {params.outdir}/metadata_target.tsv
				
				unzip target.zip -d {params.targetdir}

				shuf -n 500 {params.targetdir}/ncbi_dataset/fetch.txt > {params.targetdir}/ncbi_dataset/tmp
				mv {params.targetdir}/ncbi_dataset/tmp {params.targetdir}/ncbi_dataset/fetch.txt
				datasets rehydrate --directory {params.targetdir}
				rm target.zip
				mv {params.targetdir}/ncbi_dataset/data/GC?_*/*.fna {params.targetdir}/
				rm -rf {params.targetdir}/ncbi_dataset {params.targetdir}/README.md

				find {params.targetdir} -type f -name 'GC*.fna' -exec bash -c 'mv "$1" "$(dirname "$1")/$(basename "$1" | cut -d"_" -f1,2).fna"' _ {{}} \;
			fi
		
		elif [[ {params.download_target} == "no" ]]
		then
			echo "Target genomes provided"
			cp -r {params.user_target} {output.genomedir}
		else
			echo "No input found"
		fi
		
		ls {params.targetdir}/*.fna | awk -F "/" '{{ print $NF }}' | sed 's/.fna//g' > {output.genomelist}
		
		"""


checkpoint download_offtarget:
	input:
		infile = "myseq.fasta",
	output:
		genomedir = directory(GENOMES_OFFTARGET),
		genomelist = OUTDIR + "genomes_offtarget.txt"
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
	threads: 4
	shell:
		"""
		if [[ {params.download_offtarget} == "yes" ]]
		then
			if [[ {params.domain} == "viral" && {params.use_assembly} == "no" ]]
			then
				datasets download virus genome taxon {params.spid} --complete-only --filename offtarget.zip
				dataformat tsv virus-genome --package offtarget.zip > {params.outdir}/metadata_offtarget.tsv

				unzip offtarget.zip -d {params.offtargetdir}
				rm offtarget.zip
				seqkit split --by-id -O {params.offtargetdir} {params.offtargetdir}/ncbi_dataset/data/genomic.fna
				ls {params.offtargetdir}/*.fna | parallel -j {threads} 'rename "genomic.part_" "" "{{}}"'  
				rm -rf {params.offtargetdir}/ncbi_dataset
				rm {params.offtargetdir}/README.md
			else
				datasets download genome taxon {params.spid} \
				--assembly-source {params.db} \
				--assembly-version latest \
				--assembly-level {params.assembly_level} \
				--include genome \
				--dehydrated \
				--filename offtarget.zip
				
				unzip offtarget.zip -d {params.offtargetdir}

				shuf -n 1000 {params.offtargetdir}/ncbi_dataset/fetch.txt > {params.offtargetdir}/ncbi_dataset/tmp
				mv {params.offtargetdir}/ncbi_dataset/tmp {params.offtargetdir}/ncbi_dataset/fetch.txt
				datasets rehydrate --directory {params.offtargetdir}
				dataformat tsv genome --package offtarget.zip > {params.outdir}/metadata_offtarget.tsv
				
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
		
		ls {params.offtargetdir}/*.fna | awk -F "/" '{{ print $NF }}' | sed 's/.fna//g' > {output.genomelist}
		
		"""

rule prepare_primers:
	input:
		primers = PRIMERS,
	output:
		primers = OUTDIR + "primers.txt",
		primers_expand = OUTDIR + "primers_expand.txt",
		probes = OUTDIR + "probes.fasta",
		status = OUTDIR + "status/prepare_primers.txt",
	params:
		outdir = OUTDIR,
	threads: 8
	conda: "../envs/download.yaml"
	shell:
		"""
		mkdir -p {params.outdir}
		awk -F "\\t" '{{ print $1, toupper($2), toupper($3), toupper($4) }}' OFS="\\t" {input.primers} > {params.outdir}/primerstmp.txt
		python scripts/expand_iupac.py {params.outdir}/primerstmp.txt {output.primers_expand}
		
		cut -f2-4 {output.primers_expand} > {output.primers}
		cut -f2,5 {output.primers_expand} | sed 's/^/>/1' | sed 's/\\t/\\n/g' > {output.probes}
		rm {params.outdir}/primerstmp.txt

		touch {output.status}
		"""