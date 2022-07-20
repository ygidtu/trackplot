## Docker Usage

### Command line

```bash
docker run --rm -v $PWD:$PWD --user $(id -u):$(id -g) ygidtu/sashimi --help
```

- `-v`, `--volumn`: mount the working directory to docker container
- `--user`: prevent docker read and write file using root privileges
- the rest usage please check [Command line usage](./command_line_usage.md)

