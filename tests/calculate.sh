#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find -name * | grep -v _dis | xargs md5sum 
