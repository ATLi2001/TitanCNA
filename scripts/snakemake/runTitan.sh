#!/bin/bash

# get the reference and liftover files from bucket
# mkdir /data/
# echo "Make /data/ directory done!"
# gsutil cp -r gs://bucket1q2w/jpliu/titancna/tcga/refgrch38/cytoBand_hg38.txt /data/
# echo "Copy ref data done!"
# gsutil cp -r gs://bucket1q2w/jpliu/titancna/tcga/liftover /data/
# echo "Copy liftover data done!"
####
## Mount a disk on existing container direcotry will miss all the files in container.
### Here I copy them to the mounted disk for process.
cp -r /home/usr/* /mnt/data
cd /mnt/data/titancna/scripts/snakemake/
snakemake -s TitanCNA.snakefile -j 100 -p --config bamMeta=$bamMeta \
samplesList=$samplesList \
pairList=$pairList \
run=$run \
token=$token

cp -r results/titan/hmm/* $OUTPUT_DIR
cp -r logs/* $OUTPUT_DIR
cp -r data/gdc/ $OUTPUT_DIR


