#!/usr/bin/env bash

if [ $# -lt 1 ]; then
  echo -e "\nUsage: getcontainers.sh [containers directory]\n"
else
  containers_dir="$1"

  for container in $(grep -w 'container' configs/containers/base.config | sed 's/=/\t/g' | cut -f2); do
    singularity \
      pull ${containers_dir}$(echo ${container/docker:\/\//} | sed 's/[:\/]/-/g' | sed 's/"//g' | sed "s|'||g").img \
      $(echo ${container} | sed 's/"//g' | sed "s|'||g");
  done
fi
