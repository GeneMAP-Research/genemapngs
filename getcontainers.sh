#!/usr/bin/env bash

containers_dir="/scratch4/awonkam1/containers/"

for container in $(grep -w 'container' configs/containers/base.config | sed 's/=/\t/g' | cut -f2); do
  singularity \
    pull ${containers_dir}$(echo ${container/docker:\/\//} | sed 's/[:\/]/-/g' | sed 's/"//g' | sed "s|'||g").img \
    $(echo ${container} | sed 's/"//g' | sed "s|'||g");
done
