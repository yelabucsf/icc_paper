#!/bin/bash

bcftools query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
    dat/processed/cm_il1b/reqtl/vcfs/chrALL.dose.flt.bcf \
    > dat/processed/cm_il1b/reqtl/vcfs/chrALL.info.vcf

bed_fn="./dat/processed/cm_il1b/reqtl/vcfs/chrALL.info.bed"
echo -e "#chrom\tchromStart\tchromEnd\tname\tscore\tstrand" > $bed_fn
bcftools query \
    -f '%CHROM\t%POS0\t%POS\t%ID\t1000\t.\n' \
    dat/processed/cm_il1b/reqtl/vcfs/chrALL.dose.flt.bcf \
    >> $bed_fn
