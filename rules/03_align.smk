checkpoint split_amplicons_target:
	input:
		infile = rules.collate_results.output.status,
	output:
		alndir = directory(OUTDIR + "alignments_target"),
		status = OUTDIR + "status/split_target.txt",
	conda: "../envs/seqkit.yaml"
	threads: 8
	params:
		outdir = OUTDIR,
		target = OUTDIR + "ispcr_target/target_amplicons.fasta",	
	shell:
		"""
		if [[ -s {params.target} ]]
		then
			seqkit split --by-id -j {threads} --two-pass --update-faidx \
				--id-regexp '^\S+\s+(\S+)' --out-dir {output.alndir} {params.target}
			ls {output.alndir}/*.fasta | parallel -j {threads} 'rename "target_amplicons.part_" "" "{{}}"'
			ls {output.alndir}/*.fasta | awk -F "/" '{{ print $NF }}' | sed 's/.fasta//g' > {params.outdir}/target_primerlist.txt
		else
			mkdir {output.alndir}
		fi
		touch {output.status}
		"""

checkpoint split_amplicons_offtarget:
	input:
		infile = rules.collate_results.output.status,
	output:
		alndir = directory(OUTDIR + "alignments_offtarget"),
		status = OUTDIR + "status/split_offtarget.txt",
	conda: "../envs/seqkit.yaml"
	threads: 8
	params:
		outdir = OUTDIR,
		offtarget = OUTDIR + "ispcr_offtarget/offtarget_amplicons.fasta",		
	shell:
		"""
		if [[ -s {params.offtarget} ]]
		then
			seqkit split --by-id -j {threads} --two-pass --update-faidx \
				--id-regexp '^\S+\s+(\S+)' --out-dir {output.alndir} {params.offtarget}
			ls {output.alndir}/*.fasta | parallel -j {threads} 'rename "offtarget_amplicons.part_" "" "{{}}"'
			ls {output.alndir}/*.fasta | awk -F "/" '{{ print $NF }}' | sed 's/.fasta//g' > {params.outdir}/offtarget_primerlist.txt
		else
			mkdir {output.alndir}
		fi
		touch {output.status}
		"""


rule check_split:
	input:
		target = get_primersets_t,
		offtarget = get_primersets_ot,
	output:
		status = OUTDIR + "status/check_split.txt",
	shell:
		"""
		touch {output.status}
		"""


rule align_amplicons_target:
	input:
		rules.check_split.output.status,
		#checkpoints.split_amplicons_target.get(**wildcards).output[2],
	output:
		#alignment_t = directory(OUTDIR + "alignments_target"),
		status = OUTDIR + "status/align_amplicon_target.txt",
	conda: "../envs/align.yaml"
	threads: 16
	params:
		alndir = OUTDIR + "alignments_target",
		ispcrdir = OUTDIR + "ispcr_target",
		primerlist = OUTDIR + "target_primerlist.txt",
	shell:
		"""
		if [[ -s {params.primerlist} ]]
		then
			cat {params.primerlist} | \
				parallel -j {threads} \
				'ginsi --thread 4 \
				{params.alndir}/"{{}}".fasta > {params.alndir}/"{{}}".aln'
		fi
		touch {output.status}
		"""


rule align_amplicons_offtarget:
	input:
		rules.check_split.output.status,
		#checkpoints.split_amplicons_offtarget.get(**wildcards).output[2],
	output:
		#alignment_ot = directory(OUTDIR + "alignments_offtarget"),
		status = OUTDIR + "status/align_amplicon_offtarget.txt",
	conda: "../envs/align.yaml"
	threads: 16
	params:
		alndir = OUTDIR + "alignments_offtarget",
		ispcrdir = OUTDIR + "ispcr_offtarget",
		primerlist = OUTDIR + "offtarget_primerlist.txt",
	shell:
		"""
		if [[ -s {params.primerlist} ]]
		then
			cat {params.primerlist} | \
				parallel -j {threads} \
				'ginsi --thread 4 \
				{params.alndir}/"{{}}".fasta > {params.alndir}/"{{}}".aln'
		fi
		touch {output.status}
		"""
		