# Single-cell-miRNA-prediction
## Introduction
Here we present the computation methods to infer the single-cell miRNA expression from the scRNA expression.
Two methods are showed here: Rank-based method and Regression method.<br>
<br>
The rank-based method used the miRNA database-TargetScan and cancer data from TCGA to find the RNA targets. Then the miRNA expression was calculated based on the rank
of its target across all the genes in one cell. The regression method directly uses the miRNA-RNA information learning from the bulk RNA and miRNA data by the machine learning model. In our analysis, the rank-based method could get more robust result with different quality of scRNA data while the regression method could achieve the best accuracy when the input scRNA data is imputed.<br>
<br>
The code in the git hub is based on the analysis on breast data and TNBC data in the TCGA database.

