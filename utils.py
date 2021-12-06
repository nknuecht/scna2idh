# ==============================================================================
__author__ = "Nicholas Nuechterlein"
__license__ = "MIT"
__maintainer__ = "Nicholas Nuechterlein"
# ==============================================================================
import pandas as pd

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
        p_value_dict[str((class1, class2))] = results.p_value

    p_value_df = pd.DataFrame.from_dict(p_value_dict, orient='index').sort_values(0)
    return p_value_df
