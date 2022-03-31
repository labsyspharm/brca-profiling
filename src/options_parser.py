import os
from optparse import OptionParser
import subprocess

brca_dir = os.path.dirname(os.path.dirname(__file__))
data_dir = os.path.join(brca_dir, "data")
results_dir = os.path.join(os.getcwd(), "brca_results")
if os.path.isdir(results_dir):
    pass
else:
    subprocess.run(["mkdir", results_dir])

# Option parser
###############
def arguments():
    parser = OptionParser()
    parser.add_option("-b", "--baseline-file",  dest="baseline_data", type="str", default="%s/rnaseq_log2rpkm.csv"%data_dir,
                    help="baseline data - RNAseq, mass-spectometry, phospho mass-spectometry, kinase activity")
    parser.add_option("-r", "--response-file", dest="response_data", type="str", default="%s/grmetrics.csv"%data_dir,
                    help="dose-response for drug-cell combinations, columns headers - 'GR_AOC', 'sigma_GR_AOC', 'GRmax', 'sigma_GR_AOC', 'GEC', 'sigma_GEC50'")
    parser.add_option("-g", "--gene-list", dest="gene_list", type="str", default="%s/gene_list.txt"%data_dir,
                    help="subset of genes on which the predictions are made")
    parser.add_option("-c", "--cell-list", dest="cell_list", type="str", default="%s/cell_list.txt"%data_dir,
                    help="subset of genes on which the predictions are made")
    parser.add_option("-d", "--drug", dest="drug", type="str", default="abemaciclib", help="drug name")
    parser.add_option("-m", "--metric", dest="metric", type="str", default="GR_AOC",
                    help="metric to be predicted, possible values: GR_AOC, GRmax, GEC50")
    parser.add_option("-p", "--parameters", dest="params", type="str", default="%s/randomforest_params.txt"%data_dir,
                    help="file containing the parameters for RandomForestRegressor")
    parser.add_option("-t", "--prediction_type", dest="prediction_type", type="str", default="predict_genes",
                    help="switch between evaluating feature importance, estimating model accuracy, or recurssive feature elimination - predict_genes, estimate_accuracy, recurssive_elimination, estimate_error")
    parser.add_option("-n", "--num_pairs", dest="complexity", type="int", default=5,
                    help="Computational complexity, if None then AUC is evaluated over all possible cell pairs, by default each cell is paired with 5 other cells.")
    parser.add_option("-f", "--num_features", dest="num_features", type="int", default=300,
                    help="Number of features for feature pruning")
    parser.add_option("-s", "--split", dest="test_set", type="str", default="%s/test_lines.csv"%data_dir,
                    help="file with train and test lines for recurssive feature elimination")
    parser.add_option("-o", "--output", dest="output", type="str", default=results_dir,
                    help="output folder")
    options, args = parser.parse_args()

    return options, args
