#! /bin/bash

blast_db="/usr/share/ncbi/data/swissprot"
blast_evalue='1e-10'
query_file=$1

#echo $query_file

#echo "${blastp_bin} -db ${blast_db} -evalue '${blast_evalue}' -query ${query_file}"

blastp -db ${blast_db} -evalue "${blast_evalue}" -query ${query_file}
