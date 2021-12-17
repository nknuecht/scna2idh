# ==============================================================================
__author__ = "Nicholas Nuechterlein"
__license__ = "MIT"
__maintainer__ = "Nicholas Nuechterlein"
# ==============================================================================
import pandas as pd
import numpy as np

from lifelines.statistics import multivariate_logrank_test

def get_pairwise_p_values(kmf_df):
    '''
    Return a pandas dataframe containing the pairwise logrank
    pvalues of group survival curves
    '''

    vals = np.unique(kmf_df[kmf_df.columns[1]])
    com_list = []

    for i in range(len(vals)):
        for j in range(len(vals)):
            if i != j and (vals[j], vals[i]) not in com_list:
                com_list.append((vals[i], vals[j]))

    p_value_dict = {}
    for class1, class2 in com_list:
        pair_idxs = kmf_df.loc[kmf_df[kmf_df.columns[1]].isin([class1, class2])].index
        kmf_temp_df = kmf_df.loc[pair_idxs]
        col0, col1, col2 = kmf_temp_df.columns[0], kmf_temp_df.columns[1], kmf_temp_df.columns[2]
        results = multivariate_logrank_test(kmf_temp_df[col0],  kmf_temp_df[col1], kmf_temp_df[col2])
        p_value_dict[str((class1, class2))] = {'p-value':results.p_value}

    p_value_df = pd.DataFrame.from_dict(p_value_dict, orient='index').sort_values('p-value')
    return p_value_df

## Nicholas Nuechterlein
import numpy as np
import pandas as pd
from sklearn import metrics

def get_metrics_helper(preds, probs, labels):
    '''
    return AUC, balanced accuracy, F1, precision, recall, kappa score, and
    mathews correlation coeffient for binary classifer predictions
    '''
    # auc score
    if len(np.unique(labels)) == 2:
        auc_score = metrics.roc_auc_score(np.asarray(labels), probs)
    else:
        auc_score = np.nan

    # average accuracy
    class0_idxs = labels.loc[labels == 0].index.tolist()
    class1_idxs = labels.loc[labels == 1].index.tolist()
    class0_acc = np.mean(labels.loc[class0_idxs] == preds.loc[class0_idxs])
    class1_acc = np.mean(labels.loc[class1_idxs] == preds.loc[class1_idxs])
    avg_acc = np.mean([class0_acc, class1_acc])

    # f1 score
    f1 = metrics.f1_score(np.asarray(labels), preds)

    # precision and recall
    precision = metrics.precision_score(np.asarray(labels), preds)
    recall = metrics.recall_score(np.asarray(labels), preds)

    # kappa and
    kappa = metrics.cohen_kappa_score(y1=np.asarray(labels), y2=preds)
    MCC = metrics.matthews_corrcoef(y_true=np.asarray(labels), y_pred=preds)

    # confusion matrix
    if len(np.unique(labels)) == 2:
        cm = metrics.confusion_matrix(np.asarray(labels), preds)
        cm_df = pd.DataFrame(data=cm, columns=['Pred 0', 'Pred 1'], index=['True 0', 'True 1'])
    else:
        cm_df = np.nan

    # store and return results in a dictonary
    metric_dict = {'auc':auc_score, 'avg_acc':avg_acc, 'f1':f1,  'precision':precision, 'recall':recall,
                   'kappa':kappa, 'MCC':MCC}
    return metric_dict, cm_df

def get_metrics(preds, probs, labels, thresh=0.95):
    '''
    return metrics of binary classifier for high, low, and all model prediction confidence levels
    '''
    metric_all_dict, cm_all_df = get_metrics_helper(preds=preds, probs=probs, labels=labels)


    high_probs_idxs = probs.loc[(probs < (1 - thresh)) | (probs > thresh)].index
    low_probs_idxs = probs.loc[(probs >= (1 - thresh)) & (probs <= thresh)].index

    metric_high_dict, cm_high_df = get_metrics_helper(preds=preds.loc[high_probs_idxs],
                                                      probs=probs.loc[high_probs_idxs],
                                                      labels=labels.loc[high_probs_idxs])
    metric_low_dict, cm_low_df = get_metrics_helper(preds=preds.loc[low_probs_idxs],
                                                    probs=probs.loc[low_probs_idxs],
                                                    labels=labels.loc[low_probs_idxs])

    metric_dict = {'All samples':metric_all_dict,
                   'High confidence samples':metric_high_dict,
                   'Low confidence samples':metric_low_dict}
    cm_dict = {'All samples':cm_all_df,
               'High confidence samples':cm_high_df,
               'Low confidence samples':cm_low_df}

    return metric_dict, cm_dict
