#!/bin/bash
. VERSION && echo "tags: $TAGS"
for t in $TAGS; do
   docker image tag images.canfar.net/uvickbos/find_moving:latest images.canfar.net/uvickbos/find_moving:$t
done
unset TAGS
docker image list images.canfar.net/uvickbos/find_moving
