#!/bin/bash

awk -F "\t" 'BEGIN {{
				comp["A"] = "T"; comp["T"] = "A"; comp["C"] = "G"; comp["G"] = "C"
			}}
			$2 != 4 && $2 != 20 {{
				if ($2 ~ /16/) {{
					# If the alignment is on the reverse strand (FLAG 16), reverse complement the sequence
					seq = $10
					rev_comp = ""
					for (i = length(seq); i > 0; i--) {{
						rev_comp = rev_comp comp[substr(seq, i, 1)]
					}}
					$10 = rev_comp
				}}
				# Print the line with the sequence either reversed or not
				print
			}}' OFS="\t" "$1" | sed 's/^r_//g' > "$1".rc