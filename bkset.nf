nextflow.enable.dsl=2

params.ng   = 50                    // Number of genes
params.ns   = 30                    // Number of sets
params.out  = 'background'          // Write sets to this directory
params.pfx  = 'bkset'               // Prefix to use in the output file

// Data file to use for sampling feature names
params.platform = 'rna'

process randgenes {
    container "${params.contPfx}labsyspharm/brca-profiling:${params.contVers}"
    publishDir "${params.out}", mode: 'move', saveAs: {f -> "${params.pfx}-${params.ng}-${index}.txt"}

    input: val(index)
    output: path('bkset.txt')
    
    script:
    blf = '/app/data/' + (params.platform == 'ms' ? 'mass_spec.csv' : 'rnaseq_log2rpkm.csv')
    """
    cut -d ',' -f 1 $blf | tail -n +2 | shuf | head -n ${params.ng} > bkset.txt
    """
}

workflow {
    Channel.of(1..params.ns) | randgenes
}
