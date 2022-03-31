import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from options_parser import arguments
from datasets import brca_data
from random_forest import feature_importance

# OptionParser
options, args = arguments()
# Breast cancer data
dff = brca_data(full_data=True)
train_cells, test_cells = brca_data(train_test=True)

# Random forest parameters
##########################
dfparam = pd.read_csv(options.params, sep="\t")
param_name = dfparam["parameter"].tolist()
value = dfparam["value"].tolist()
for i in range(dfparam.shape[0]): exec("%s=%d"%(param_name[i], value[i]))

def feature_prune(df):
    '''
    Recurssive feature elimination
    '''
    df = df[df.index.isin(train_cells)]
    itr=1
    features = []
    drop_itr = []
    while df.shape[1]>=300:
        dfimp = feature_importance(df)
        drop_features = dfimp.tail(int(dfimp.shape[0]*0.25)).index.tolist()
        for feat in drop_features:
            features.append(feat)
            drop_itr.append(int(itr))
        df = df.drop(columns=drop_features)
        itr+=1
    df = df.drop(columns=["%s"%options.metric,"sigma_%s"%options.metric])
    keep_features = df.columns.tolist()
    for feat in keep_features:
        features.append(feat)
        drop_itr.append(np.nan)
    dffeat = pd.DataFrame(list(zip(features,drop_itr)), columns=["features", "drop_iteration"])
    return dffeat

def train_test_eval(df, dffeat, train_cells, test_cells):
    drop_features = dffeat.dropna()["features"].tolist()
    df = df.drop(columns=drop_features)
    dftrain = df[df.index.isin(train_cells)]
    dftest = df[df.index.isin(test_cells)]
    ytrain = dftrain[options.metric].values
    ytrain_err = dftrain["sigma_%s"%options.metric].values
    ytest = dftest[options.metric].values
    ytest_err = dftest["sigma_%s"%options.metric].values
    dftrain = dftrain.drop(columns=["%s"%options.metric, "sigma_%s"%options.metric])
    Xtrain = dftrain.values
    dftest = dftest.drop(columns=["%s"%options.metric, "sigma_%s"%options.metric])
    Xtest = dftest.values

    random_forest = RandomForestRegressor(n_estimators=n_estimators, max_depth=max_depth,
                                          random_state=random_state)
    random_forest.fit(Xtrain, ytrain)
    ypred = random_forest.predict(Xtest)
    num_rp = 0
    auc = 0
    dfout = pd.DataFrame(columns=["cell_line1", "cell_line2", "measured1", "measured2", "predicted1", "predicted2"])
    for i in range(len(ytest)):
        for j in range(i+1, len(ytest)):
            if abs(ytest[i]-ytest[j])>max(ytest_err[i], ytest_err[j]):
                num_rp+=1
            else:
                continue
            if (ypred[i]-ypred[j])*(ytest[i]-ytest[j])>0: auc+=1
            df1 = pd.DataFrame({"cell_line1":[dftest.index[i]], "cell_line2":[dftest.index[j]],
                                "measured1":[ytest[i]], "measured2": [ytest[j]],
                                "predicted1":[ypred[i]], "predicted2":[ypred[j]]})
            dfout = pd.concat([dfout, df1], ignore_index=True)
    auc = float(auc/num_rp)
    return auc, dfout

dffeat = feature_prune(dff)
auc, dfout = train_test_eval(dff, dffeat, train_cells, test_cells)
dffeat.to_csv("%s/%s/%s_rfe.csv"%(options.output, options.drug, options.drug), index=False)
dfout.to_csv("%s/%s/%s_pairs.csv"%(options.output, options.drug, options.drug), index=False)
with open("%s/%s/%s_auc.txt"%(options.output, options.drug, options.drug),"w") as outFile:
    outFile.write(str(auc))
outFile.close()
