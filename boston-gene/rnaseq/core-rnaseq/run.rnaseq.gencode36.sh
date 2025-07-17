nextflow run \
    nf-core/rnaseq \
    -profile docker \
    --input /home/michael/cheng-project/boston-gene/rnaseq/core-rnaseq/cheng.rna.samplesheet.csv \
    --aligner star_rsem \
    --fasta /home/michael/references/gencode-v36/GRCh38.primary_assembly.genome.fa \
    --gtf /home/michael/references/gencode-v36/gencode.v36.annotation.gtf \
    --star_index /home/michael/references/gencode-v36/star_index/ \
    -c nextflow.medium.config \
    --outdir /home/michael/cheng-project/boston-gene/rnaseq/core-rnaseq/results-v36/ \
   --save_reference
