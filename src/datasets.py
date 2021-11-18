import pandas as pd

def load_response(gr_values=False):
    dfgr = pd.read_csv("../data/grmetrics.csv")
    dfgrvals = pd.read_csv("../data/grvalues.csv", index_col=0)
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
