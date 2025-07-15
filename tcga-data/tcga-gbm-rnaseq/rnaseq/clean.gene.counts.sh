#!/bin/bash
#awk 'NR==2 {for(i=1;i<=NF;i++) if($i=="gene_id") g=i; else if($i=="tpm_unstranded") t=i; print $g "\t" $t; next} $1 ~ /^ENSG/ {print $g "\t" $t}' input.txt > gene_tpm.tsv

for dir in */ ; do
  infile=$(find "$dir" -maxdepth 1 -name "*.rna_seq.augmented_star_gene_counts.tsv" | head -n 1)
  if [[ -f "$infile" ]]; then
    base=$(basename "$infile" .rna_seq.augmented_star_gene_counts.tsv)
    outfile="${base}.rna_seq.augmented_star_gene_counts_tpm_clean.tsv"
    awk 'NR==2 {for(i=1;i<=NF;i++) if($i=="gene_id") g=i; else if($i=="tpm_unstranded") t=i; print $g "\t" $t; next} $1 ~ /^ENSG/ {print $g "\t" $t}' "$infile" > "$outfile"
    echo "✅ Processed: $infile → $outfile"
  else
    echo "❌ Skipped: No matching file in $dir"
  fi
done
