# Single-cell-miRNA-prediction
## Introduction
Here we present the computation methods to infer the single-cell miRNA expression from the scRNA expression.
Two methods are showed here: Rank-based method and Regression method.<br>
<br>
The rank-based method used the miRNA database-TargetScan and cancer data from TCGA to find the RNA targets. Then the miRNA expression was calculated based on the rank
of its target across all the genes in one cell. The regression method directly uses the miRNA-RNA information learning from the bulk RNA and miRNA data by the machine learning model. In our analysis, the rank-based method could get more robust result with different quality of scRNA data while the regression method could achieve the best accuracy when the input scRNA data is imputed.<br>
<br>
The code in the git hub is based on the analysis on breast data and TNBC data in the TCGA database.
## Requirement
* R 4.0.2
* Seurat
* SAVER
* TCGA
* TargetScan
## Getting started
1. Prepare the requirement in your local computer
2. Download the method code you want to apply in the Github
3. Input the scRNA data into the code step by step
## Note
* The breast simulation data and the corresponding imputation data is too large to store in the Github, you can download from http://predscmiRNA.chengqi.site/Data/Simulation/Raw/data_simulation(10000).Rdata and http://predscmiRNA.chengqi.site/Data/Simulation/Imputation/data_imputed.RData.
* Currently, the code and analysis is for the breast cancer data, you could also modified the code and apply to the cancer data you are interested in.
## Online service
We also offer the [online service](http://predscmiRNA.chengqi.site/PredscmiRNA.php) to infer the scmiRNA.
## Citation
If you use this data, tool or code, please considering citing:
Zhao, C., Cheng, Q., Xie, W., Xu, J., Xu, S., Wang, Y., & Feng, W. (2022). Methods for predicting single-cell miRNA in breast cancer. Genomics, 114(3), 110353. Advance online publication. https://doi.org/10.1016/j.ygeno.2022.110353
