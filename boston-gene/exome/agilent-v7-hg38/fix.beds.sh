(base) michael@head-node:~/cheng-project/boston-gene/exome/agilent-v7-hg38$ cat S31285117_Padded.bed | grep -v '^track' | grep -v '^browser' | cut -f1-3 | awk 'OFS="\t"{print $1,$2,$3}' | less
