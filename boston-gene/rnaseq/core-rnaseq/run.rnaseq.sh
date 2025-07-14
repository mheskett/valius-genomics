nextflow run \
    nf-core/rnaseq \
    -profile docker \
    --input /home/michael/cheng-project/boston-gene/rnaseq/core-rnaseq/cheng.rna.samplesheet.csv \
    --aligner star_rsem \
    --fasta /home/michael/references/gencode/GRCh38.primary_assembly.genome.fa \
    --gtf /home/michael/references/gencode/gencode.v44.annotation.gtf \
    --star_index /home/michael/references/gencode/star-genome/ \
    -c nextflow.config \
    --outdir /home/michael/cheng-project/boston-gene/rnaseq/core-rnaseq/results/
