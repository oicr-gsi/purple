#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

#find all files, return their md5sums to std out
unzip -p OCT_011633_Ki_P_OCT_011633-TS.purple.zip OCT_011633_Ki_P_OCT_011633-TS.purple.cnv.gene.tsv | md5sum
unzip -p OCT_011633_Ki_P_OCT_011633-TS.purple.zip OCT_011633_Ki_P_OCT_011633-TS.purple.qc  | md5sum
unzip -p OCT_011633_Ki_P_OCT_011633-TS.purple.zip OCT_011633_Ki_P_OCT_011633-TS.purple.cnv.somatic.tsv |md5sum

