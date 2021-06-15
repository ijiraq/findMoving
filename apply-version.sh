#!/bin/bash
. VERSION && echo "tags: $TAGS"
for t in $TAGS; do
   docker image tag images.canfar.net/uvickbos/findMoving:latest images.canfar.net/uvickbos/findMoving:$t
done
unset TAGS
docker image list images.canfar.net/uvickbos/findMoving
