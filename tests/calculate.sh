#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

ls | sed 's/.*\.//' | sort | uniq -c
find . -name '*.tsv'  | xargs md5sum | sort
for i in $(find . -name '*.vcf.gz'); do if [ -f $i ]; then zcat $i | grep -v ^# | md5sum; fi; done
