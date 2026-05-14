
DEVNAME = find_moving

NAME = images.canfar.net/uvickbos/$(DEVNAME)
VERSION = 2.1

build: dependencies Dockerfile
	docker build -t $(NAME):$(VERSION) -f Dockerfile .

push: build
	docker push $(NAME):$(VERSION) 

dependencies: 

init:
	mkdir -p build

.PHONY: clean
clean:
	\rm -rf build
