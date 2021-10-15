# ==============================================================================
__author__ = "Nicholas Nuechterlein"
__license__ = "MIT"
__maintainer__ = "Nicholas Nuechterlein"
# File Description: Reformat example data
# ==============================================================================

import pandas as pd

# load tsv SCNA matrix of [genes x patients] and save as csv file of [patients x genes]
scna_path = '../data/TCGA.GBMLGG.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz'
print('Loading Xena Data')
scna_df = pd.read_csv(scna_path, sep='\t', index_col=0).T
print('Saving reformatted Xena data of')
print('Shape:\t', scna_df.shape)
print('Patients:\t', scna_df.index[:3])
print('Genes:\t', scna_df.columns[:3])
scna_df.to_csv('../data/Xena_GBMLGG_GISTIC_Scores.csv')
