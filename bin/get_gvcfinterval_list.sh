#!/bin/bash

####################################################
# make overlapping intervals of N base pairs for   #
# chromosomes larger than N base pairs. Otherwise, #
# print out chromosome length as only interval.    #
####################################################

if [ $# -lt 4 ]; then
  echo "Usage: $(basename $0) [vcf] [chunk-size] [overlap] [output-dir]"
else
  vcf=/data/awonkam1/kesoh/wgs/data/genemapwgs-674-hg38.PASS.vcf.gz
  chunk_size=20000000
  overlap=1000
  
  vcf=$1; chunk_size=$2; overlap=$3; out=$4
  
  bcftools \
      index \
      --stats \
      ${vcf} \
      > contigs.txt
  
  awk '
    {
      if($1 ~ /^chr/ && length($1) <= 5) {
        print $1,$2
      } else if(!($1 ~ /^chr/) && length($1) <= 2) {
        print $1,$2
      }
    }' contigs.txt | \
  awk \
    -v interval=${chunk_size} \
    -v overlap=${overlap} '
      {
        if( $2 <= interval ) {
          print $1,"0",$2
        } else{
          print $1,"0",(interval+overlap);
          for( i=interval; i<=$2; i+=interval ) {
  	  if( (i+interval) < $2 ) {
              print $1,(i-overlap),(i+interval+overlap)
            } else{
                print $1,(i-overlap),$2
            }
          }
        }
      }' \
    > .interval_list
  
  while read interval; do
    echo $interval > ${out}/$(echo ${interval} | sed 's/ /./g').bed
  done < .interval_list
  rm .interval_list
fi
