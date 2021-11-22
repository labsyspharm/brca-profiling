import os
import subprocess
import pandas as pd
import numpy as np
from optparse import OptionParser
import warnings
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
from lpocv.lpocv import LeavePairOut

# Option parser
###############
parser = OptionParser()
parser.add_option("-b", "--baseline-file",  dest="baseline_data", type="str", default="../data/rnaseq_log2rpkm.csv",
                  help="baseline data - RNAseq, mass-spectometry, phospho mass-spectometry, kinase activity")
parser.add_option("-r", "--response-file", dest="response_data", type="str", default="../data/grmetrics.csv",
                  help="dose-response for drug-cell combinations, columns headers - 'GR_AOC', 'sigma_GR_AOC', 'GRmax', 'sigma_GR_AOC', 'GEC', 'sigma_GEC50'")
parser.add_option("-g", "--gene-list", dest="gene_list", type="str", default="../data/gene_list.txt",
                  help="subset of genes on which the predictions are made")
parser.add_option("-c", "--cell-list", dest="cell_list", type="str", default="../data/cell_list.txt",
                  help="subset of genes on which the predictions are made")
parser.add_option("-d", "--drug", dest="drug", type="str", default="Abemaciclib/LY2835219", help="drug name")
parser.add_option("-m", "--metric", dest="metric", type="str", default="GR_AOC",
                  help="metric to be predicted, possible values: GR_AOC, GRmax, GEC50")
parser.add_option("-p", "--parameters", dest="params", type="str", default="../data/randomforest_params.txt",
                  help="file containing the parameters for RandomForestRegressor")
parser.add_option("-t", "--prediction_type", dest="prediction_type", type="str", default="predict_genes",
                  help="switch between evaluating feature importance and estimating model accuracy - predict_genes, estimate_accuracy")
parser.add_option("-n", "--num_pairs", dest="complexity", type="int", default=5,
                  help="Computational complexity, if None then AUC is evaluated over all possible cell pairs, by default each cell is paired with 5 other cells.")
parser.add_option("-f", "--num_features", dest="num_features", type="int", default=50,
                  help="Number of features for feature pruning")
parser.add_option("-o", "--output", dest="output", type="str", default="../results/",
                  help="output folder")
options, args = parser.parse_args()

# Random forest parameters
##########################
dfparam = pd.read_csv(options.params, sep="\t")
param_name = dfparam["parameter"].tolist()
value = dfparam["value"].tolist()
for i in range(dfparam.shape[0]): exec("%s=%d"%(param_name[i], value[i]))

# Data
########
input_cells = open(options.cell_list, "r").read().splitlines()
input_cells = [ic.replace("-", "").upper() for ic in input_cells]
dfbase = pd.read_csv(options.baseline_data, index_col=0)
base_cells = dfbase.columns.tolist()
genes = open(options.gene_list, "r").read().splitlines()
dfgr = pd.read_csv(options.response_data)
dfgr = dfgr[dfgr.agent==options.drug]
dfgr.index=dfgr.cell_line
gr_cells = dfgr.index.tolist()
cells = set(input_cells).intersection(set(base_cells).intersection(set(gr_cells)))
dfbase = dfbase[cells]
dfgr = dfgr[dfgr.index.isin(cells)]
dfc = pd.concat([dfbase[dfbase.index.isin(genes)].T,
                 dfgr[[options.metric, "sigma_%s"%options.metric]]], axis=1).dropna()
dff = pd.concat([dfbase.T,
                 dfgr[[options.metric, "sigma_%s"%options.metric]]], axis=1).dropna()

#####################################################################################################################################
def feature_importance(df):
    '''
    Function to evaluate feature importance of input feature set. 
    Parameters:
    ###########
    df: pandas dataframe, shape=(rows,columns) 
        rows: features (genes)
        columns: samples + growth rate metrics 
        must contain a column "METRIC" and its coressponding "sigma_METRIC"
    Returns:
    ########
    dfimp: pandas dataframe
        shape = (rows,columns) rows: features (genes)
                               columns: feature importance, rho, pval of Spearman rank correlation between
                               gene expression and growth rate metric across all cells
    '''
    X = df.drop(columns=[options.metric, "sigma_%s"%options.metric]).values
    y = df[options.metric].values
    random_forest = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth,
                                          random_state=random_state)
    random_forest.fit(X, y)                   
    importances = random_forest.feature_importances_
    indices = np.argsort(importances)[::-1]
    feature_labels = df.drop(columns=[options.metric, "sigma_%s"%options.metric]).columns
    dfimp = pd.DataFrame(list(zip(feature_labels[indices], importances[indices])),
                         columns=['features', options.drug])
    dfimp.index = dfimp.features.tolist()
    dfimp["spearman_rho"] = [None]*dfimp.shape[0]
    dfimp["spearman_pval"] = [None]*dfimp.shape[0]
    for feature in dfimp.index:
        rho, pval = spearmanr(df[feature], df[options.metric])
        dfimp.loc[feature, "spearman_rho"] = rho
        dfimp.loc[feature, "spearman_pval"] = pval
    return dfimp


def feature_pruning(df, num_features):
    '''
    Iterative feature pruning for data-driven feature selection, drops 25% of the least important features in every iteration
    Parameters:
    ##########
    df: pandas dataframe, shape=(rows,columns) 
        rows: features (genes)
        columns: samples + growth rate metrics 
        must contain a column "METRIC" and its coressponding "sigma_METRIC"
 
    num_features: int, default=50
            number of features to be returned
    '''
    dfimp = feature_importance(df)
    while dfimp.shape[0] > (options.num_features+1):
        drop_features = dfimp.tail(int(0.25*dfimp.shape[0])).index.tolist()
        df = df.drop(columns=drop_features)
        dfimp = feature_importance(df)
    return dfimp


def compute_auc(df):
    '''
    Random forest regressor.
    Parameters:
    ##########
    df: pandas dataframe, shape=(rows,columns) 
        rows: features (genes)
        columns: samples + growth rate metrics 
        must contain a column "METRIC" and its coressponding "sigma_METRIC"
    Returns:
    ########
    auc: float
        Estimate of model accuracy
    dfout: pandas dataframe
        Measured and predicted values in each cross-validation fold.
    '''
    X = df.drop(columns=[options.metric, "sigma_%s"%options.metric]).values
    y = df[options.metric].values
    cell_line = dfc.index.tolist()
    yerr = df["sigma_%s"%options.metric].values
    auc=0
    itr=0
    dfout = pd.DataFrame()
    num_pairs=options.complexity
    for train, test in LeavePairOut().split(X, y, 2.0*yerr, num_pairs=num_pairs):
        random_forest = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth,
                                              random_state=random_state)
        random_forest.fit(X[train], y[train])
        ypred = random_forest.predict(X[test])
        itr+=1
        df0=pd.DataFrame(data=[["cv_iteration%s"%itr, cell_line[0], y[test[0]], ypred[0]]],
                         columns=["CV-fold", "cell_line", "measured", "predicted"])
        df1=pd.DataFrame(data=[["cv_iteration%s"%itr, cell_line[1], y[test[1]], ypred[1]]],
                         columns=["CV-fold", "cell_line", "measured", "predicted"])

        dfout = dfout.append(df0, ignore_index=True)
        dfout = dfout.append(df1, ignore_index=True)
        if (ypred[0]-ypred[1])*(y[test[0]]-y[test[1]]) > 0:
            auc+=1
    auc = float(auc/itr)
    return auc, dfout

####################################################################################################################################

if options.prediction_type=="predict_genes":
    print ("Predicting drivers for %s (%s cell lines)"%(options.drug, dfc.shape[0]))
    dfimp = feature_importance(dfc)
    dfimp.to_csv("%s%s_imp.csv"%(options.output, options.drug.split("/")[0]), index=False)

elif options.prediction_type=="feature_pruning":
    print("Predicting drivers (data driven) for %s (%s cell lines)"%(options.drug, dff.shape[0]))
    dfimp = feature_pruning(dff, options.num_features)
    dfimp.to_csv("%s%s_ddfs.csv"%(options.output, options.drug.split("/")[0]), index=False)

elif options.prediction_type=="estimate_accuracy":
    print ("Estimating accuracy of Random Forest model for %s (%s cell lines)"%(options.drug, dfc.shape[0]))
    auc, dfout = compute_auc(dfc)
    dfout.to_csv("%s%s.csv"%(options.output, options.drug.split("/")[0]), index=False)
    with open("%s%s.txt"%(options.output, options.drug.split("/")[0]),"w") as outFile:
        outFile.write(str(auc))
    outFile.close()
