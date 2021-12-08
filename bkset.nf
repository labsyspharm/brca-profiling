nextflow.enable.dsl=2

params.ng   = 50                    // Number of genes
params.ns   = 30                    // Number of sets
params.out  = 'background'          // Write sets to this directory
params.pfx  = 'bkset'               // Prefix to use in the output file

// Data file to use for sampling feature names
params.baselineFile = 'rnaseq_log2rpkm.csv'

process randgenes {
    container "${params.contPfx}labsyspharm/brca-profiling:${params.contVers}"
    publishDir "${params.out}", mode: 'move', saveAs: {f -> "${params.pfx}-${params.ng}-${index}.txt"}

    input: val(index)
    output: path('bkset.txt')
    
    """
    cut -d ',' -f 1 /app/data/${params.baselineFile} | \
      tail -n +2 | shuf | head -n ${params.ng} > bkset.txt
    """
}

workflow {
    Channel.of(1..params.ns) | randgenes
}
