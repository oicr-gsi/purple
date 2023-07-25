#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

ls | sed 's/.*\.//' | sort | uniq -c
find . -name '*.tsv'  | xargs md5sum | sort
find . -name '*.vcf.gz' | xargs zcat | grep -v ^# | md5sum | sort
