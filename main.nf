nextflow.enable.dsl=2

// Expected params
// .in - directory that contains gene signatures, with each signature
//       in a separate .txt file, listing one gene name per line
params.out          = 'results.csv'
params.baselineFile = 'rnaseq_log2rpkm.csv'

process accuracy {
    container "${params.contPfx}labsyspharm/brca-profiling:${params.contVers}"
    
    input:
    tuple path(genes), val(drug)
    file(cl)   // cell lines

    output:
    tuple val("${genes.getBaseName()}"), path('*.txt')

    script:
    carg = params.containsKey('cellList') ? "-c $cl" : '-c /app/data/cell_list.txt'
    """
    python /app/src/random_forest.py -t estimate_accuracy \
      -b /app/data/${params.baselineFile} \
      $carg -d '$drug' -g $genes -o ./
    """
}

workflow {
    cell_lines = params.containsKey('cellList') ? file(params.cellList) : ''
    sigs  = Channel.fromPath("${params.in}/*.txt")
    drugs = Channel.of(params.drugs).flatten()
    inputs = sigs.combine(drugs)

    f = file("${params.out}")
    if( f.exists() ) error "File ${params.out} already exists"
    f << "Signature,Drug,AUC\n"
    
    accuracy(inputs, cell_lines)
        .map{sig, f -> "$sig,${f.getBaseName().split('_auc').head()},${f.text}"}
        .subscribe{ f << "$it\n" }
}
