nextflow.enable.dsl=2

process accuracy {
    container 'labsyspharm/brca-profiling:1.0.0'
    
    input:
    val(drug)
    file(cl)   // cell lines

    output:
    path('*.txt')

    script:
    carg = params.containsKey('cellLines') ? "-c $cl" : '-c /app/data/cell_list.txt'
    """
    python /app/src/random_forest.py -t estimate_accuracy \
      $carg -d $drug -o ./ \
      -b /app/data/rnaseq_log2rpkm.csv \
      -r /app/data/grmetrics.csv \
      -g /app/data/gene_list.txt \
      -p /app/data/randomforest_params.txt
    """
}

process aggregate {
    input:
    val(aucs)

    output:
    path('test.csv')

    script:
    """
    echo Drug,AUC > test.csv
    echo "$aucs" >> test.csv
    """
}

workflow {
    drugs = Channel.of(params.drugs).flatten()
    cell_lines = params.containsKey('cellLines') ?
        file(params.cellLines) : ''
        
    accuracy(drugs, cell_lines)
        .map{it -> "${it.getBaseName()},${it.text}"}
        .reduce{a, b -> "$a\n$b"} |
        aggregate
}
