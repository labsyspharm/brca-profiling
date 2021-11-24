# brca-profiling

**Prerequisites:** [Docker](https://docs.docker.com/get-docker/), [nextflow](https://www.nextflow.io/)

Place all your signature `*.txt` files into a single directory (e.g., `my_signatures`), then run the following:

```
nextflow run labsyspharm/brca-profiling --in my_signatures
```

The script will generate a `results/` directory and populate it with AUC values from evaluating each signature against each drug.
The output directory can be controlled with `--out`, e.g.,

```
nextflow run labsyspharm/brca-profiling --in my_signatures --out /path/to/results
```

Signatures can be evaluated on a subset of cell lines. Make a `cl.txt` file that lists which cell lines should be considered (one per line), then feed it to the script with `--cell-lines`:

```
nextflow run labsyspharm/brca-profiling --in my_signatures --cell-lines cl.txt
```
