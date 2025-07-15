awk 'BEGIN {OFS=FS="\t"} {sub(/\.[0-9]+$/, "", $1); print}' $1 > output
