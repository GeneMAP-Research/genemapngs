####################################################
# make overlapping intervals of N base pairs for   #
# chromosomes larger than N base pairs. Otherwise, #
# print out chromosome length as only interval.    #
####################################################

awk '
  {
    if(\$1 ~ /^chr/ && length(\$1) <= 5) {
      print \$1,\$2
    } else if(!(\$1 ~ /^chr/) && length(\$1) <= 2) {
      print \$1,\$2
    }
  }' ${contigs} | \
awk \
  -v interval=${params.chunk_size} \
  -v overlap=${params.overlap} '
    {
      if( \$2 <= interval ) {
        print \$1,"0",\$2
      } else{
        print \$1,"0",(interval+overlap);
        for( i=interval; i<=\$2; i+=interval ) {
	  if( (i+interval) < \$2 ) {
            print \$1,(i-overlap),(i+interval+overlap)
          } else{
              print \$1,(i-overlap),\$2
          }
        }
      }
    }' \
  > .interval_list


while read interval; do
  echo \$interval > \$(echo \${interval} | sed 's/ /./1; s/ /_/g').bed
done < .interval_list

