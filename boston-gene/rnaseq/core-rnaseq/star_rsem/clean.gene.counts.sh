####
awk '{print $1 "\t" $NF}' rsem.merged.gene_tpm.tsv > rsem.merged.gene_tpm.cleaned.tsv

