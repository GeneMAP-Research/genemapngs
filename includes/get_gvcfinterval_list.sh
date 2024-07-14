#!/usr/bin/env bash

######################################################
# make non-overlapping intervals of 5,000,000 bp for #
# chromosomes larger than 5,000,000. Otherwise, just #
# print out chromosome length as only interval.      #
# HLA contigs will be automatically excluded. To do  #
# HLA typing, use the SNP2HLA or equivalent tool     #
######################################################

if [ $# -lt 3 ]; then
  echo -e "Usage: get_interval_list.sh [gvcf-file] [interval-size] [output-prefix]"
else
  gvcf=$1  #"/scratch/eshkev001/projects/wes/higenes/careni/hg38/gvcfs/OS.P001_92401_S25_val.bqsr.g.vcf.gz"
  intval=$2
  out=$3
  
  zgrep '##contig' ${gvcf} | \
      sed 's/[=,>]/\t/g' | \
      cut -f3,5 | \
      grep -v '^HLA' | \
      awk \
        -v intvl=${intval} '
          { 
	    if( $2 <= intvl ) { 
	      print $1,"0",$2 
            } else{ 
	        for( i=0; i<=$2; i+=intvl ) { 
	          if( i+(intvl-1) < $2 ) { 
		    print $1,i,i+(intvl-1)
	          } else{ 
		      print $1,i,$2 
	          } 
	        } 
	    } 
	  }' \
  > ${out}_interval_list.txt
fi

#while read interval; do
#    echo $interval > $(echo ${interval} | sed 's/[:*]/_/g' | sed 's/ /_/g').bed
#done < .interval_list

