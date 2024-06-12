#!/usr/bin/env bash

work="$(pwd)/workdir/"

./genemapngs.sh test

nextflow -c test.config run workflow/test.nf -w ${work}
