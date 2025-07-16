#!/bin/bash

for dir in */ ; do
  infile=$(find "$dir" -maxdepth 1 -name "*.rna_seq.augmented_star_gene_counts.tsv" | head -n 1)
  if [[ -f "$infile" ]]; then
    base=$(basename "$infile" .rna_seq.augmented_star_gene_counts.tsv)
    outfile="${base}.rna_seq.augmented_star_gene_counts_tpm_clean.tsv"
    
    awk 'NR==2 {
        for(i=1;i<=NF;i++) {
            if($i=="gene_id") g=i;
            else if($i=="tpm_unstranded") t=i;
        }
        print "gene_id\ttpm_unstranded";
        next
    }
    $1 ~ /^ENSG/ {
        sub(/\.[0-9]+$/, "", $g);  # remove .version from ENSG ID
        print $g "\t" $t
    }' "$infile" > "$outfile"
    
    echo "✅ Processed: $infile → $outfile"
  else
    echo "❌ Skipped: No matching file in $dir"
  fi
done
