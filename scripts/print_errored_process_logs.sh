#!/usr/bin/env bash

grep 'terminated with an error exit status' .nextflow.log | grep 'Execution is retried' | sed -E 's|^.+? - \[([a-z0-9]+/[a-z0-9]+)\] .+$|work/\1\*|g' | xargs -I @ bash -c 'echo =================================================== @/.command.sh ===================================================; cat @/.command.sh; echo =================================================== @/.command.log ===================================================; cat @/.command.log;echo;echo;echo;echo'
