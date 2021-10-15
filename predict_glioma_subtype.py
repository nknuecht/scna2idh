# ==============================================================================
__author__ = "Nicholas Nuechterlein"
__license__ = "MIT"
__maintainer__ = "Nicholas Nuechterlein"
# File Description: This file contains the code for predicting adult glioma molecular subtype
# ==============================================================================

from argparse import ArgumentParser

import pandas as pd
import numpy as np

from pycaret.classification import load_model, predict_model

def predict(args):

    scna_path = args.scna_path
    gene_loc_path = args.gene_loc_path
    model_path = args.model_path
    outfile = args.outfile
    threshold = args.threshold

    # load scna data and gene location data
    print('Loading SCNA data . .  . ')
    scna_df = pd.read_csv(scna_path, index_col=0)
    gene_loc_df = pd.read_csv(gene_loc_path, index_col=0)

    # exclude sex chromosomes and small chromosome arms (21p, 22p)
    chr_arms = ['1p', '1q', '2p', '2q', '3p', '3q', '4p', '4q', '5p', '5q', '6p', '6q', '7p', '7q',
                '8p', '8q', '9p', '9q', '10p', '10q', '11p', '11q', '12p', '12q', '13q', '14q', '15q',
                '16p', '16q', '17p', '17q', '18p', '18q', '19p', '19q', '20p', '20q', '21q', '22q']
    gene_loc_df = gene_loc_df.loc[gene_loc_df['chr_arm'].isin(chr_arms)]

    # find genes with chromosome arm annotations in the gene_loc file
    common_genes = [x for x in scna_df.columns if x in gene_loc_df.index]

    # Treat homozyous deletions as single-copy loss & high-level amplifications as low-level amplifications
    scna_df = scna_df[common_genes].replace({-2:-1, 2:1})

    # average over chromosome arms
    chrarm_scna_df = pd.concat([gene_loc_df['chr_arm'], scna_df.T], axis=1, join='inner').groupby('chr_arm').mean().T

    # create a dataframe that only considers losses for the 1p/19q-codeletion screen
    scna_loss_df = scna_df[common_genes].replace({-2:-1, 2:1})

    # average over chromosome arms
    chrarm_scna_loss_df = pd.concat([gene_loc_df['chr_arm'], scna_df.replace({1:0}).T],
                                    axis=1, join='inner').groupby('chr_arm').mean().T

    ## Stage 1 ##
    # screen for 1p/19q-codeletion (oligodendroglioma prediction)
    print('Predicting 1p/19q-codeletions . .  . ')
    thresh = 0.85
    oligo_pred_idxs = chrarm_scna_loss_df.loc[(chrarm_scna_loss_df['1p'] < -thresh)
                                              & (chrarm_scna_loss_df['1q'] >= -thresh)
                                              & (chrarm_scna_loss_df['19p'] >= -thresh)
                                              & (chrarm_scna_loss_df['19q'] < -thresh)].index.tolist()
    astro_pred_idxs = [x for x in chrarm_scna_loss_df.index if x not in oligo_pred_idxs]

    # make oligo prediction df
    oligo_pred_df = pd.DataFrame(data=np.ones(len(oligo_pred_idxs))*2, index=oligo_pred_idxs, columns=['Pred'])
    oligo_pred_df['Score'] = 1
    oligo_pred_df = oligo_pred_df[['Score', 'Pred']]

    ## Stage 2 ##
    # load classifier
    idh_classifier = load_model(model_path)
    # SCNA data from tumors without predicted 1p/19q-codeletions
    astro_chrarm_scna_df = chrarm_scna_df.drop(index=oligo_pred_idxs)

    # predict IDH-mutations
    print('Predicting IDH-mutations for samples without predicted 1p/19q-codeletions . .  . ')
    astro_pred_df = predict_model(idh_classifier, data = astro_chrarm_scna_df)
    astro_pred_df = astro_pred_df[['Score', 'Label']].rename(columns={'Label':'Pred'})

    ## merge predictions ##
    prediction_df = oligo_pred_df.append(astro_pred_df)
    prediction_df['Pred'] = prediction_df['Pred'].replace({'0.0':'IDH-Wildtype Glioma', '1.0':'IDH-Mutant Astroctyoma',
                                                           2:'Oligodendroglioma'})

    # assign low confidence predictions with a "Uncertain" label
    print('Retaining predictions with confidence over'+str(threshold)+' . .  . ')
    low_confidence_preds_idxs = prediction_df.loc[prediction_df['Score' < float(threshold)]].index
    prediction_df.loc[low_confidence_preds_idxs, 'Pred'] = 'Uncertain'

    if outfile:
        prediction_df.to_csv(outfile)

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--scna_path', '-s',default='./data/TCGA.GBMLGG.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz')
    parser.add_argument('--gene_loc_path', '-g', default="./data/gistic_cytoband_chr_arm_23109x4.csv")
    parser.add_argument('--model_path', '-m', default='./models/model_10-14-2021')
    parser.add_argument('--outfile', '-o', default="./data/predictions.csv")
    parser.add_argument('--threshold', '-t', type=float, default=0.7)
    args = parser.parse_args()

    predict(args=args)
