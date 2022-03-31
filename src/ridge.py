import pandas as pd
from sklearn.linear_model import Ridge
from options_parser import arguments
from lpocv.lpocv import LeavePairOut
from datasets import brca_data

# OptionParser
options, args = arguments()

# Breast cancer data
dfc = brca_data()

#AUC ridge regression
def auc_ridge(df):
    '''
    Parameters:
    ##########
    df: pandas dataframe, shape=(rows,column"""s)
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
    for train, test in LeavePairOut().split(X, y, 2.0*yerr, num_pairs=options.complexity):
        ridge = Ridge()
        ridge.fit(X[train], y[train])
        ypred = ridge.predict(X[test])
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

auc, dfout = auc_ridge(dfc)
dfout.to_csv("%s/%s_ridge.csv"%(options.output, options.drug), index=False)
with open("%s/%s_ridge.txt"%(options.output, options.drug),"w") as outFile:
    outFile.write(str(auc))
outFile.close()
