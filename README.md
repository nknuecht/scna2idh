# scna2idh

![Predictive System](./figures/predictive_system.png)

This repository gives a system that classifies adult diffuse glioma as 1p/19q-codeleted oligodendroglioma, IDH-mutant astrocytoma, or IDH-wildtype glioblastoma using somatic copy number alteration (SCNA) data alone. This system was developed using cross-validation on 786 adult diffuse gliomas in The Cancer Genome Atlas (TCGA)  projects TCGA-LGG and TCGA-GBM, validated using three independant dataset and a holdout portion of the TCGA, and depolyed on subjects in the REMBRANDT study for whom SCNA data, but not IDH mutational status, is available. 

## Results

### Cross-validation

![cv_results](./figures/cv_results.png)

Calling 1p/19q-codeletions based a 85% threshold for loss of chromosome arms 1p, 1q, 19p, and 19q our system classifies oligodendrogliomas in all three versions of our TCGA training set (UCSC hg19, GDC hg19, GDC hg38) with 100% accuacy. Our system's second stage--an IDH-mutation classifier trained on adult astrocytic glioma--also performed well  (AUC= 0.990 +/- 0.001, MCC=0.935 +/- 0.006). We calabrated this classifier and rejected samples whose prediction was given with 70% or less confidence. This raised our model's performance (AUC= 0.992 +/- 0.001, MCC = 0.970 +/- 0.004) at the cost of only rejecting 5% of samples.

### Validation on three independent datasets and hold TCGA validation set

![val_results](./figures/val_results.png)

As shown above, an 85% threshold was optimal for the TCGA training set and optimal or nearly optimal for two independent validation sets published by Capper et al. (MCC=0.97) and Jonsson et al. (MCC=0.97), respectively. Our IDH mutation classifier performed well across three independent validation datasets and our holdout TCGA cohort of histological grade 4 patients with surrogate IDH labels not found from IDH sequencing. When evaluated on patients with model confidence greater than 70%, our IDH mutation classifier archived AUC scores greater than 0.95 on each dataset.

### REMBRANDT study predictions

![rembrandt_results](./figures/rembrandt_results.png)

Our REMBRANDT 1p/19q-codeleted oligodendroglioma screen captured the densest area of tumors with significant gene losses on chromosome arms 1p and 19q, consistent with the 1p/19q-codeletion screens on our training and validation sets. Our IDH mutation predictions on histological astrocytomas generated a dramatic survival difference between predicted IDH-wildtype gliomas and IDH-mutant astrocytomas (HR=1.51, p = 0.003), log-rank) and appeared to correctly identify histological lower-grade IDH-wildtype gliomas (OS=1.1 years) now considered to be IDH-wildtype glioblastomas. Histological glioblastomas predicted to be IDH-wildtype glioblastomas have the same median overall survival (1.1 years), and their survival trajectory was significantly worse than those of histological grade 4 tumors that were predicted IDH-mutant astrocytomas (HR=0.78, p=0.004, log-rank).

## Useage 

1. Generate copy number segmentation files. 
2. Run GISTIC 2.0.
   - Create GISTIC 2.0 input files as specificed by GenePattern's [documentation](https://www.genepattern.org/modules/docs/GISTIC_2.0).
   - GISTIC 2.0 pameters used in this study are listed in the forthcoming paper.
3. Make predictions
   - python file: 
   - jupyter-notebook: 
