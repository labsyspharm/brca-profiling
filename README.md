# brca-profiling: Evaluate the capacity of gene sets to predict drug response in breast cancer cell lines

**Prerequisites:** [Docker](https://docs.docker.com/get-docker/), [nextflow](https://www.nextflow.io/)

**Step 1**

Place all your gene sets into a single directory. Each gene set must be in its own `.txt` file, listing one gene name per line. For example,

```
mkdir genesets

echo APAF1 BAD BAK1 BAX BCL2 BCL2A1 BCL2L1 BCL2L11 BID BIK CASP3 CASP7 CASP8 \
CASP9 HRK PMAIP1 PUMA BMF BCL2L2 CYCS MCL1 | tr ' ' '\n' > genesets/apoptosis.txt

echo AKT1 AKT2 AKT3 CRKL IRS2 MTOR PIK3CA PIK3CG PIK3R1 PIK3R2 PTEN RICTOR \
RNF43 RPTOR TSC1 TSC2 | tr ' ' '\n' > genesets/pi3k.txt
```

**Step 2**

Pull the latest version of the code and run it with the directory composed in Step 1:

```
nextflow pull labsyspharm/brca-profiling
nextflow run labsyspharm/brca-profiling --in genesets
```

The script will generate a `results/` directory and populate it with AUC values from evaluating each signature against each drug.
The output directory can be controlled with `--out`, e.g.,

```
nextflow run labsyspharm/brca-profiling --in my_signatures --out /path/to/results
```

Signatures can be evaluated on a subset of cell lines. Make a `cl.txt` file that lists which cell lines should be considered (one per line), then feed it to the script with `--cell-lines`:

```
echo 184A1 AU565 BT20 BT474 BT549 CAL120 | tr ' ' '\n' > cl.txt
nextflow run labsyspharm/brca-profiling --in my_signatures --cell-lines cl.txt
```

## Running on O2

Docker is not required to run the script on O2, where Singularity will be used instead. All the necessary parameters are encapsulated inside the configuration profile. Simply specify `-profile O2` when executing the script:

```
module load java
nextflow run labsyspharm/brca-profiling --in genesets -profile O2
```
