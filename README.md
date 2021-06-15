## building 
Update VERSION
```
docker build -t images.canfar.net/uvickbos/find_moving:latest -f Dockerfile .
docker run --interactive --tty --rm  --env DISPLAY=host.docker.internal:0 --volume ${DATA_DIR}:/DATA images.canfar.net/uvickbos/find_moving:latest bash
./apply-version.sh
```

