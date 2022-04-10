nextflow.enable.dsl=2

// Expected params
// .in - directory that contains gene signatures, with each signature
//       in a separate .txt file, listing one gene name per line
params.out      = 'results.csv'
params.platform = 'rna'

process accuracy {
    container "${params.contPfx}labsyspharm/brca-profiling:${params.contVers}"
    
    input:
    tuple path(genes), val(drug)
    file(cl)   // cell lines

    output:
    tuple val("${genes.getBaseName()}"), path('*.txt')

    script:
    carg = params.containsKey('cellList') ? "-c $cl" : '-c /app/data/cell_list.txt'
    blf = '/app/data/' + (params.platform == 'ms' ? 'mass_spec.csv' : 'rnaseq_log2rpkm.csv')
    """
    python /app/src/random_forest.py -t estimate_accuracy \
      -b $blf $carg -d '$drug' -g $genes -o ./
    """
}

workflow {
    cell_lines = params.containsKey('cellList') ? file(params.cellList) : ''
    sigs  = Channel.fromPath("${params.in}/*.txt")
    drugs = Channel.of(params.drugs).flatten()
    inputs = sigs.combine(drugs)

    f = file("${params.out}")
    if( f.exists() ) error "File ${params.out} already exists"
    f << "Signature,Drug,Method,AUC\n"
    
    accuracy(inputs, cell_lines)
        .map{sig, f -> tokens = f.getBaseName().split('_');
          "$sig,${tokens[0]},${tokens[1]},${f.text}"}
        .subscribe{ f << "$it\n" }
}
