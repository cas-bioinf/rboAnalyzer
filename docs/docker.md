# Run rboAnalyzer using Docker

Docker image containing the rboAnalyzer and essential data.

## Install Docker

Install compatible docker version from https://www.docker.com/.


## rboAnalyzer container structure
- The container is at https://hub.docker.com/repository/docker/schwarzmarek/rboanalyzer
- The working directory of the rboAnalyzer in the container is `/data`. You are expected to mount your data there.
- The Rfam database is in the container and it's version is 14.1.
- The container is ready to use - no other preparation should be necessary

## Prepare for running the examples
- Pull rboAnalyzer image
    ```
    docker pull schwarzmarek/rboanalyzer:0.1
    ```

- Obtain the `examples` directory

    Obtain a copy of `examples` directory. (Download and unpack https://github.com/cas-bioinf/rboAnalyzer/releases/download/v0.1.0/example.zip)

- In the terminal, go into the unpacked `examples` directory and run docker command(s) from there.

### Run interactive docker session 
Inside the session you can run rboAnalyzer commands. Take care to point the output files to the mounted directory (otherwise they will not be saved).

The commnand below expects that you are inside the directory with the data you want to analyze (e.g. the `examples` directory) and does following:
- Runs the `/bin/bash` in the specified container (`schwarzmarek/rboanalyzer:0.1`) in the interactive mode (`-it`)
- Mounts the current directory (`source="$(pwd)"`) to `/data` directory (`target=/data`) inside the container. This will cause that any changes in the mounted directory inside the container will be propagated to the directory on host. (for other options see https://docs.docker.com/storage/)

```
docker container run -it --mount type=bind,source="$(pwd)",target=/data schwarzmarek/rboanalyzer:0.1 /bin/bash
```

Now you should see the content of the `examples` dictionary. If so, you can execute the commands from the [Example](../readme.md#Example) section.