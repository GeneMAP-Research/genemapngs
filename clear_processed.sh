#!/usr/bin/env bash

##########################################################################
# Watch a directory for new simbolic links and follow the paths to the   #
# parent (intermediate) files and delete them. Then replace the symbolic #
# links in the target directory with the actual files                    #
#                                                                        #
# Useful when running nextflow to clear up space for jobs that have      #
# completed successfully                                                 #
##########################################################################

if [ ! $# -eq 1 ]; then
  echo "Usage: clear_processed.sh [target directory]"
else
  dir=$1
  if [ ! -d $dir ]; then
    echo "Have you provided the correct target directory? Please check."
    exit 1
  fi
  #dir=/scratch/eshkev001/projects/wgs/camscd/hg38/cram/
  readlink ${dir}/* > target_dir_status.txt
  mkdir -p ${dir}/temp
  
  while [ -s target_dir_status.txt ]; do
  
    echo "Symbolic links found"
  
    for i in $(readlink ${dir}/*); do readlink $(dirname $i)/*; done > rm_processed_0.txt
    
    for m in $(for l in $(for k in $(for j in $(for i in $(readlink ${dir}/*); do readlink $(dirname $i)/*; done); do readlink $(dirname $j)/*; done); do readlink $(dirname $k)/*; done); do readlink $(dirname $l)/*; done); do readlink $(dirname $m)/*; done > rm_processed_1.txt
    
    for l in $(for k in $(for j in $(for i in $(readlink ${dir}/*); do readlink $(dirname $i)/*; done); do readlink $(dirname $j)/*; done); do readlink $(dirname $k)/*; done); do readlink $(dirname $l)/*; done > rm_processed_2.txt
    
    for k in $(for j in $(for i in $(readlink ${dir}/*); do readlink $(dirname $i)/*; done); do readlink $(dirname $j)/*; done); do readlink $(dirname $k)/*; done > rm_processed_3.txt
    
    for j in $(for i in $(readlink ${dir}/*); do readlink $(dirname $i)/*; done); do readlink $(dirname $j)/*; done > rm_processed_4.txt 
    
    cat rm_processed_{0..4}.txt | sort | uniq > rm_processed.txt
  
    for i in $(cat rm_processed.txt); do [ -e ${i} ] && rm -rf $i; done
  
    mv $(readlink ${dir}/*) ${dir}/temp

    for n in ${dir}/*; do 
        [ ! -e $i ] && rm $n; 
    done

    mv ${dir}/temp/* ${dir}/

    rmdir ${dir}/temp

    rm rm_processed*
  
    readlink ${dir}/* > target_dir_status.txt
  
  done
fi
