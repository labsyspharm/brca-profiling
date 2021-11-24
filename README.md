# brca-profiling

**Prerequisites:** [Docker](https://docs.docker.com/get-docker/), [nextflow](https://www.nextflow.io/)

Place all your signature `*.txt` files into a single directory (e.g., `my_signatures`), then run the following:

```
nextflow run labsyspharm/brca-profiling --in `my_signatures`
```

The script will generate a `results/` directory and populate it with AUC values from evaluating each signature against each drug.
