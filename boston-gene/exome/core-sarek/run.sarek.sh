# note - chatGPT gets all of this wrong. do manually
## remove ASCAT until I get the exome files established.
## might need to do custom ref file prep for cnvkit as well

nextflow run nf-core/sarek \
  -profile docker \
  --wes \
  --input /home/michael/cheng-project/boston-gene/exome/core-sarek/cheng.samplesheet.csv \
  --genome GATK.GRCh38 \
  --intervals /home/michael/cheng-project/boston-gene/exome/agilent-v7-hg38/S31285117_Padded_cleaned.bed \
  --tools haplotypecaller,mutect2,manta,cnvkit,snpeff,vep \
  --outdir /home/michael/cheng-project/boston-gene/exome/core-sarek/results/ \
  --save_mapped \
  --save_output_as_bam \
  -c /home/michael/cheng-project/boston-gene/exome/core-sarek/nextflow.config
#  -resume # only if resuming
