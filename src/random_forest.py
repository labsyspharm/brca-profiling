import os
import subprocess
import warnings
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import spearmanr
from warnings import filterwarnings
from lpocv.lpocv import LeavePairOut
from options_parser import arguments
from datasets import brca_data


# OptionParser
options, args = arguments()

# Breast cancer data
dfc = brca_data()

# Random forest parameters
##########################
dfparam = pd.read_csv(options.params, sep="\t")
param_name = dfparam["parameter"].tolist()
value = dfparam["value"].tolist()
for i in range(dfparam.shape[0]): exec("%s=%d"%(param_name[i], value[i]))

# AUC random forest regressor
def auc_random_forest(df):
    '''
    Parameters:
    ##########
    df: pandas dataframe, shape=(rows, columns)
        rows: features (genes)
        columns: samples + growth rate metrics
        must contain a column METRIC (e.g. "GR_AOC", "GR_max", "GR50", etc.) and its coressponding "sigma_METRIC"
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
    for train, test in LeavePairOut().split(X, y, 2.0*yerr, num_pairs=options.complexity):
        random_forest = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth, random_state=random_state)
        random_forest.fit(X[train], y[train])
        ypred = random_forest.predict(X[test])
        itr+=1
        df0=pd.DataFrame(data=[["cv_iteration%s"%itr, cell_line[0], y[test[0]], ypred[0]]],
                         columns=["CV-fold", "cell_line", "measured", "predicted"])
        df1=pd.DataFrame(data=[["cv_iteration%s"%itr, cell_line[1], y[test[1]], ypred[1]]],
                         columns=["CV-fold", "cell_line", "measured", "predicted"])

        dfout = pd.concat([dfout, df0], ignore_index=True)
        dfout = pd.concat([dfout, df1], ignore_index=True)
        if (ypred[0]-ypred[1])*(y[test[0]]-y[test[1]]) > 0:
            auc+=1
    print("Evaluated %s pairs using leave-pair-out cross-validation."%itr)
    auc = float(auc/itr)
    return auc, dfout

def feature_importance(df):
    '''
    Function to evaluate feature importance of input feature set. 
    Parameters:
    ###########
    df: pandas dataframe, shape=(rows,columns) 
        rows: features (genes)
        columns: samples + growth rate metrics 
        must contain a column METRIC (e.g. "GR_AOC", "GR_max", "GR50", etc.) and its coressponding "sigma_METRIC"
    Returns:
    ########
    dfimp: pandas dataframe
        shape = (rows,columns) rows: features (genes)
        columns: feature importance, rho, pval of Spearman rank correlation between gene expression and growth rate metric across all cells
    '''
    filterwarnings("ignore")
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

#####################################################

if options.prediction_type=="predict_genes":
    print ("Predicting drivers for %s (%s cell lines)."%(options.drug, dfc.shape[0]))
    dfimp = feature_importance(dfc)
    dfimp.to_csv("%s/%s_imp.csv"%(options.output, options.drug), index=False)

elif options.prediction_type=="estimate_accuracy":
    print ("Estimating accuracy of Random Forest model for %s (%s cell lines)"%(options.drug, dfc.shape[0]))
    auc, dfout = auc_random_forest(dfc)
    dfout.to_csv("%s/%s_rfr.csv"%(options.output, options.drug), index=False)
    with open("%s/%s_rfr.txt"%(options.output, options.drug),"w") as outFile:
        outFile.write(str(auc))
    outFile.close()
else:
    print("Invalid prediction type")
