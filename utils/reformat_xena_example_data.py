# ==============================================================================
__author__ = "Nicholas Nuechterlein"
__license__ = "MIT"
__maintainer__ = "Nicholas Nuechterlein"
# File Description: Reformat example data
# ==============================================================================

import pandas as pd
import os

# load tsv SCNA matrix of [genes x patients] and
scna_path = './data/TCGA.GBMLGG.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz'
print('Loading Xena Data')
print(' -> from', os.path.abspath(scna_path))
scna_df = pd.read_csv(scna_path, sep='\t', index_col=0).T

# save SCNA matrix as csv file of [patients x genes]
outfile = './data/Xena_GBMLGG_GISTIC_Scores.csv'
print('Saving reformatted Xena data of')
print(' -> to', os.path.abspath(outfile))
print('Shape:\t', scna_df.shape)
print('Patients:\t', scna_df.index[:3])
print('Genes:\t', scna_df.columns[:3])
scna_df.to_csv(outfile)
