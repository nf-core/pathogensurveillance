#!/usr/bin/env bash

# Watch the .command.log files in active docker containers

watch 'docker stats --no-stream | tail -n +2 | cut -d " " -f 1 | xargs -I {} bash -c "echo -e "---------------{}---------------";docker exec {} tail -n 7 .command.log"'
