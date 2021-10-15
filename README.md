# scna2idh

This repository gives a system that classifies adult diffuse glioma as 1p/19q-codeleted oligodendroglioma, IDH-mutant astrocytoma, or IDH-wildtype glioblastoma using somatic copy number alteration (SCNA) data alone. This system was developed using cross-validation on 786 adult diffuse gliomas in The Cancer Genome Atlas (TCGA)  projects TCGA-LGG and TCGA-GBM, validated using three independant dataset and a holdout portion of the TCGA, and depolyed on subjects in the REMBRANDT study for whom SCNA data, but not IDH mutational status, is available. 

## Useage 

1. The first step is to generate copy number segmentation files. Segmentation files can be derived from a variety of data modalities, including SNP array, DNA methylation array, targeted sequecing, whole genome sequencing, and low-pass whole genome sequecning data. 

2. Given copy number segmenation files, GISTIC 2.0 [[1]](#1) can generate gene-level SCNA calls. The GISTIC 2.0 pameters we used in this study are listed in the forthcoming paper. GenePattern's GISTIC 2.0 [documentation](https://www.genepattern.org/modules/docs/GISTIC_2.0) is a good resource for understanding these parameters and for creating input files for GISITC 2.0. Instructions for GISTIC 2.0 installation can be found [here](http://portals.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=216&p=t).

3. A mapping between genes and chromosome arms is necessary for this. The default 

4. Requirements for our predictive system

   - The jupyter-notebook `classification_system_example.ipynb` gives a working example of 1p/19q-codeletion and IDH-mutant prediction for all samples in the TCGA-LGG and TCGA-GBM datasets. 

   - A python file can be run directly from the command line as

     ```bash
     python predict_glioma_subtype.py -s {scna_path} -g {gene_loc_path} -o {outfile}
     ```

     | Name                  | Flag | Description                                                  |
     | --------------------- | ---- | ------------------------------------------------------------ |
     | scna_path             | -s   | Filepath to a csv file whose indices are sample names, columns are genes, and its values are thresholded GISTIC 2.0 scores. |
     | gene_loc_path         | -g   | Filepath to a csv file whose indices are genes. This csv file much have a column labeled "chr_arm" that lists the chromosome arm associated with each gene. |
     | outfile               | -o   | Filepath to which patient predictions, and their prediction confidence, are saved in a .csv format. |
     | threshold (optional)  | -t   | Model only returns predictions whose prediction confidence is above the specified threshold. Default threshold is 0.7. Threshold at or below 0.5 return all predictions. |
     | model_path (optional) | -m   | Path to trained model. Default path is set to latest trained model. |

     The SCNA dataframe must contain genes from all of the following chromosome arms: 1p, 1q, 2p, 2q, 3p, 3q, 4p, 4q, 5p, 5q, 6p, 6q, 7p, 7q, 8p, 8q, 9p, 9q, 10p, 10q, 11p, 11q, 12p, 12q, 13q, 14q, 15q, 16p, 16q, 17p, 17q, 18p, 18q, 19p, 19q, 20p, 20q, 21q, 22q.

     Example

      1. Download Xena TCGA glioma SCNA data

         ```bash
         wget -P ./data/ https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.GBMLGG.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
         ```

     2. Reformat Xena TCGA glioma SCNA data

        ```bash
        python utils/reformat_xena_example_data.py 
        ```

     3. Run

        ```bash
        python predict_glioma_subtype.py -s ./data/Xena_GBMLGG_GISTIC_Scores.csv -g ./data/gistic_cytoband_chr_arm_23109x4.csv -o ./data/tcga_preds.csv 
        ```

     Note: in our study, we used TCGA SCNA data we computed from copy number segmenation files downloaded from the GDC, not data preprocessed GISTIC scores hosted by UCSC's Xena. We use Xena data here for ease. 



## References
<a id="1">[1]</a> Mermel CH, Schumacher SE, Hill B, Meyerson ML, Beroukhim R, Getz G (2011) GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biol 12: R41 Doi 10.1186/gb-2011-12-4-r41

