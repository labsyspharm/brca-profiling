import pandas as pd
from options_parser import arguments

options, args = arguments()

def brca_data(train_test=False, full_data=False):
    # input cells over which the predsictive model is built
    input_cells = open(options.cell_list, "r").read().splitlines()
    input_cells = [ic.replace("-", "").upper() for ic in input_cells]
    print("Parsed %s user-provided cell line names"%len(input_cells))

    # baseline data
    dfbase = pd.read_csv(options.baseline_data, index_col=0)
    base_cells = dfbase.columns.tolist()
    msg = "Loaded %s containing %s features and %s cell lines"
    print(msg%(options.baseline_data, dfbase.shape[0], dfbase.shape[1]))

    # list of genes (or proteins) which would serve as the features for predictive modeling
    genes = open(options.gene_list, "r").read().splitlines()
    print("Parsed %s user-provided gene names"%len(genes))

    # GR-metrics data
    dfgr = pd.read_csv(options.response_data)
    dfgr = dfgr[dfgr.generic_name==options.drug]
    dfgr.index=dfgr.cell_line
    gr_cells = dfgr.index.tolist()
    print("Loaded GR values for %s cell lines"%len(gr_cells))

    # cells that are present in both, baseline as well GR-data
    cells = list(set(input_cells).intersection(set(base_cells).intersection(set(gr_cells))))
    dfbase = dfbase[cells]
    dfgr = dfgr[dfgr.index.isin(cells)]
    print("%s cell lines are in common to all sources"%len(cells))

    # combine the baseline and dose response data into one dataFrame
    dfc = pd.concat([dfbase[dfbase.index.isin(genes)].T,
                    dfgr[[options.metric, "sigma_%s"%options.metric]]], axis=1).dropna()
    dff = pd.concat([dfbase.T,
                    dfgr[[options.metric, "sigma_%s"%options.metric]]], axis=1).dropna()

    # train-test sets, with roughly 80-20 split, for rescursive feature elimination
    dftest = pd.read_csv(options.test_set, index_col=0)
    test_cells = dftest[dftest["test"]==1].index.tolist()
    train_cells = list(set(cells)-set(test_cells))
    if train_test is False:
        if full_data is False:
            return dfc
        else:
            return dff
    else:
        return train_cells, test_cells

def load_response(gr_values=False):
    dfgr = pd.read_csv("../data/grmetrics.csv")
    dfgrvals = pd.read_csv("../data/grvalues.csv")
    if gr_values is False: 
        return dfgr
    else:
        return dfgrvals

def load_baseline():
    dfrna = pd.read_csv("../data/rnaseq_log2rpkm.csv", index_col=0)
    dfms = pd.read_csv("../data/mass_spec.csv", index_col=0)
    dfpms = pd.read_csv("../data/phospho_phase1.csv", index_col=0)
    return dfrna, dfms, dfpms

def load_meta():
    dfcells = pd.read_csv("../data/cell_lines_metadata.csv")
    dfdrugs = pd.read_csv("../data/agents_metadata.csv")
    return dfcells, dfdrugs
