# brca-profiling: Evaluate the capacity of gene sets to predict drug response in breast cancer cell lines

**Prerequisites:** [Docker](https://docs.docker.com/get-docker/), [nextflow](https://www.nextflow.io/)

**Step 1: Define gene sets**

Place all your gene sets into a single directory. Each gene set must be in its own `.txt` file, listing one gene name per line. For example,

```
mkdir genesets

echo APAF1 BAD BAK1 BAX BCL2 BCL2A1 BCL2L1 BCL2L11 BID BIK CASP3 CASP7 CASP8 \
CASP9 HRK PMAIP1 PUMA BMF BCL2L2 CYCS MCL1 | tr ' ' '\n' > genesets/apoptosis.txt

echo AKT1 AKT2 AKT3 CRKL IRS2 MTOR PIK3CA PIK3CG PIK3R1 PIK3R2 PTEN RICTOR \
RNF43 RPTOR TSC1 TSC2 | tr ' ' '\n' > genesets/pi3k.txt
```

**Step 2: Run the container on all gene sets**

Pull the latest version of the code and run it with the directory composed in Step 1:

```
nextflow pull labsyspharm/brca-profiling
nextflow run labsyspharm/brca-profiling --in genesets
```

The script will evaluate all signatures in the input directory against all drugs and write the resulting AUC values to `results.csv` by default.
The output filename can be controlled with `--out`, e.g.,

```
nextflow run labsyspharm/brca-profiling --in genesets --out my_output.csv
```

Signatures can be evaluated on a subset of cell lines. Make a `cl.txt` file that lists which cell lines should be considered (one per line), then feed it to the script with `--cell-list`:

```
echo 184A1 AU565 BT20 BT474 BT549 CAL120 | tr ' ' '\n' > cl.txt
nextflow run labsyspharm/brca-profiling --in genesets --cell-list cl.txt
```

By default, the script evaluates gene sets against RNAseq data. To evaluate against Mass Spec (MS) and phopsho-MS, use `--baseline-file` to point the script to another data file (note that the filenames are resolved w.r.t. `data/` in this repository):
```
nextflow run labsyspharm/brca-profiling --in genesets --baseline-file mass_spec.csv
nextflow run labsyspharm/brca-profiling --in genesets --baseline-file phospho_phase1.csv
```

## Evaluating against background sets

The repository includes a script for generating signatures that comprise randomly-selected genes. By default, the script generates 30 sets of 50 genes each and outputs them to `background/` subdirectory with a prefix `bkset`. The user can overwrite these values using `--ns`, `--ng`, `--out`, and `--pfx` respectively. For example,

```
nextflow run labsyspharm/brca-profiling/bkset.nf --ns 10 --ng 20 --out mysets --pfx bg
```

will generate 10 sets of 20 randomly-selected genes each and write these to `mysets/`, prefixing each filename with `bg`. The directory can then be immediately provided to the main script for evaluation:

```
nextflow run labsyspharm/brca-profiling --in mysets
```

As with gene set evaluation, background generation samples from the rnaseq space of gene names by default. This behavior can be overwritten with `--baseline-file` as above:
```
nextflow run labsyspharm/brca-profiling/bkset.nf --baseline-file mass_spec.csv --pfx ms-bg
nextflow run labsyspharm/brca-profiling/bkset.nf --baseline-file phospho_phase1.csv --pfx pms-bg
```

## Running on O2

Docker is not required to run the script on O2, where Singularity will be used instead. All the necessary parameters are encapsulated inside the configuration profile. Simply specify `-profile O2` when executing the script:

```
module load java
nextflow run labsyspharm/brca-profiling --in genesets -profile O2
```

