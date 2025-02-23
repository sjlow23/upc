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
			
		cut -f1,3 {params.outdir}/probes_expand.txt | sort | uniq | 
			awk 'BEGIN {{ prev=0 ; count=1 }} {{ if (prev==$1) count++; else {{ count=1;;prev=$1 }} print $1, $2, "probe"count }}' OFS="\\t" | 
			awk -F "\\t" '{{ print ">"$1"_"$3"\\n"$2 }}' > {output.probes}
			
		rm {params.outdir}/probestmp.txt
	
		touch {output.status}
		"""