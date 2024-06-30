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
		spid = TARGET_SP_TAXID if "TARGET_SP_TAXID" in globals() else "None",
		taxid = TARGET_TAXID if "TARGET_TAXID" in globals() else "None",
		targetdir = GENOMES_TARGET,
		db = DB,
		assembly_level = ASSEMBLY_LEVEL,
		domain = DOMAIN,
		outdir = OUTDIR,
	threads: 8
	shell:
		"""		
		if [[ {params.download_target} == "yes" ]]
		then
			if [[ {params.spid} != "None" ]]
			then
				ncbi-genome-download \
				--section {params.db} \
				--flat-output \
				--formats fasta \
				--output-folder {output.genomedir} \
				--metadata-table {params.outdir}/metadata_target.tsv \
				--assembly-levels {params.assembly_level} \
				--species-taxids {params.spid} \
				--parallel {threads} \
				--retries 5 \
				{params.domain}
			elif [[ {params.taxid} != "None" ]]
			then
				ncbi-genome-download \
				--section {params.db} \
				--flat-output \
				--formats fasta \
				--output-folder {output.genomedir} \
				--metadata-table {params.outdir}/metadata_target.tsv \
				--assembly-levels {params.assembly_level} \
				--taxids {params.taxid} \
				--parallel {threads} \
				--retries 5 \
				{params.domain}
			else
				echo "No target genomes found"
			fi

			find {params.targetdir} -type f -name "*.gz" | parallel -j {threads} gunzip

			find {params.targetdir} -type f -name "*.fna" > {params.outdir}/target_ori.txt
			find {params.targetdir} -type f -name "*.fna" | awk -F "/" '{{ print $NF }}' | \
				awk -F "_" '{{ print $1"_"$2".fna" }}' > {params.outdir}/target_rename.txt
			cat {params.outdir}/target_rename.txt | awk '{{ print "{params.targetdir}/"$0"" }}' > {params.outdir}/target_new.txt

			paste {params.outdir}/target_ori.txt {params.outdir}/target_new.txt | xargs -L1 mv

			rm {params.outdir}/target_ori.txt {params.outdir}/target_rename.txt {params.outdir}/target_new.txt

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
		spid = OFFTARGET_SP_TAXID if "OFFTARGET_SP_TAXID" in globals() else "None",
		taxid = OFFTARGET_TAXID if "OFFTARGET_TAXID" in globals() else "None",
		offtargetdir = GENOMES_OFFTARGET,
		db = DB,
		assembly_level = ASSEMBLY_LEVEL,
		domain = DOMAIN,
		outdir = OUTDIR,
	threads: 8
	shell:
		"""		
		if [[ {params.download_offtarget} == "yes" ]]
		then
			if [[ {params.spid} != "None" ]]
			then
				ncbi-genome-download \
				--section {params.db} \
				--flat-output \
				--formats fasta \
				--output-folder {output.genomedir} \
				--metadata-table {params.outdir}/metadata_offtarget.tsv \
				--assembly-levels {params.assembly_level} \
				--species-taxids {params.spid} \
				--parallel {threads} \
				--retries 5 \
				{params.domain}
			elif [[ {params.taxid} != "None" ]]
			then
				ncbi-genome-download \
				--section {params.db} \
				--flat-output \
				--formats fasta \
				--output-folder {output.genomedir} \
				--metadata-table {params.outdir}/metadata_offtarget.tsv \
				--assembly-levels {params.assembly_level} \
				--taxids {params.taxid} \
				--parallel {threads} \
				--retries 5 \
				{params.domain}
			else
				echo "No off-target genomes found"
				exit 1
			fi
		
			find {params.offtargetdir} -type f -name "*.gz" | parallel -j {threads} gunzip

			find {params.offtargetdir} -type f -name "*.fna" > {params.outdir}/offtarget_ori.txt
			find {params.offtargetdir} -type f -name "*.fna" | \
				awk -F "/" '{{ print $NF }}' | awk -F "_" '{{ print $1"_"$2".fna" }}' > {params.outdir}/offtarget_rename.txt
			cat {params.outdir}/offtarget_rename.txt | awk '{{ print "{params.offtargetdir}/"$0"" }}' > {params.outdir}/offtarget_new.txt

			paste {params.outdir}/offtarget_ori.txt {params.outdir}/offtarget_new.txt | xargs -L1 mv

			rm {params.outdir}/offtarget_ori.txt {params.outdir}/offtarget_rename.txt {params.outdir}/offtarget_new.txt

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
		probes = OUTDIR + "probes.fasta",
		status = OUTDIR + "status/prepare_primers.txt",
	threads: 8
	shell:
		"""
		cut -f1-3 {input.primers} > {output.primers}
		cut -f1,4 {input.primers} | sed 's/^/>/1' | sed 's/\\t/\\n/g' > {output.probes}
		touch {output.status}
		"""