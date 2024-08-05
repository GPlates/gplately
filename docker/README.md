### 👉 pull the docker image from Docker Hub or ghcr.io

- `docker pull ghcr.io/gplates/gplately`
- `docker pull gplates/gplately`

###  👉 build the docker image

- step 1: `git clone --depth 1 --branch master https://github.com/GPlates/gplately.git`
- step 2: `cd gplately`
- step 3: `docker build -f docker/Dockerfile -t gplately .`

Note: if errors occur, try build with **--no-cache** option

### 👉 run the docker container

```docker run --rm -ti -p 8888:8888 -v `pwd`:/workspace/my_stuff gplately```

### 👉 build and push docker images to Docker Hub and ghcr.io

- option 1: create and push a tag (will trigger the github action)
- option 2: 
    - step 1: create a branch `build-docker-image` (based on master or whatever branch you like)
    - step 2: update the version info in docker/build-docker-image-version.txt
    - step 3: push the change in step 2 (will trigger the github action)
    - note: this will also create the `latest` tag

