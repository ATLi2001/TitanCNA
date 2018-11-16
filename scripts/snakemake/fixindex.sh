#!/bin/bash

myarray=(`find $1 -maxdepth 1 -name "*.bam.bai"`)
if [ ${#myarray[@]} -gt 0 ]; then 
    echo "bam.bai exist" 
else 
    find $1 -maxdepth 1 -name "*.bai" -exec \
    bash -c 'mv $1 ${1%.bai}.bam.bai' - '{}' \;
fi