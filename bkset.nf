nextflow.enable.dsl=2

params.ng  = 50           // Number of genes
params.ns  = 30           // Number of sets
params.out = 'background' // Write sets to this directory

process randgenes {
    container "${params.contPfx}labsyspharm/brca-profiling:${params.contVers}"
    publishDir "${params.out}", mode: 'move', saveAs: {f -> "bkset-${params.ng}-${index}.txt"}

    input: val(index)
    output: path('bkset.txt')
    
    """
    cut -d ',' -f 1 /app/data/rnaseq_log2rpkm.csv | \
      tail -n +2 | shuf | head -n ${params.ng} > bkset.txt
    """
}

workflow {
    Channel.of(1..params.ns) | randgenes
}
