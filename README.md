# OSCAA
A Two-Dimensional Gaussian Mixture Model for Copy Number Variation Association Analysis

## Author
Xuanxuan Yu, Xizhi Luo, Guoshuai Cai, Feifei Xiao

## Description
The One-Stage CNV-disease Association Analysis (OSCAA) method utilizes a two-dimensional GMM model to assess the association between CNV and disease risk in a known CNVR with a one-stage framework. The two-dimensional GMM model is factorized into a signal model, a phenotype model, and a copy number model. To better differentiate samples with different copy numbers, we model the joint distribution of the first two PCs for samples with the same copy number in the signal model. The phenotype model comprises a generalized linear model (GLM) or a logistic model fitted for continuous or binary phenotype measurements, respectively. In addition, the phenotype model offers three hypotheses: the hypothesis of linear relationship between disease risk and copy number, the hypothesis of whether samples with abnormal copy number states are more/less likely to develop a disease, and the hypothesis of whether sample with deletions of copy numbers are more/less likely to develop a disease. The copy number model estimates the proportion of each copy number state and enables the evaluation of the potential impact of covariates on the proportions. Although we demonstrate the model based on SNP array data, our method can be extended to sequencing data with appropriate data normalization procedures.


## Install OSCAA
```
library(devtools)
install_github("....")
```

# Data
```
library("OSCAA")
data(exmaple.data)
dim(example.data$example.lrr)
[1] 2829   15
example.data$position
start end
1     15

example.data$example.lrr[1:6,1:6]
                                17853        80996      6145       6146       733221     6147
X4787234217_R01C01.Log.R.Ratio  0.182890500  0.1213354  0.2504343  0.1104460  0.1308962  0.1118831
X4787234217_R02C01.Log.R.Ratio -0.299842900 -0.3343200 -0.3939640 -0.5059705 -0.2845776 -0.1981067
X4787234217_R03C01.Log.R.Ratio  0.009729712 -0.1223645  0.1303647  0.1416307 -0.1090477  0.1928931
X4784244200_R01C01.Log.R.Ratio -0.423113000 -0.4317191 -0.5583948 -0.7634997 -0.5352160 -0.1343689
X4784244200_R02C01.Log.R.Ratio -0.249242800 -0.4356031 -0.2522257 -0.2990728 -0.1104633 -0.1861957
X4784244200_R03C01.Log.R.Ratio -0.211350600 -0.5990200 -0.4808382 -0.7339458 -0.3549831 -0.1585388

head(example.data$clinic.info[,c("age","Gender.char")])
     age                Gender.char
2671 50.85753           M
2410 52.87671           M
1400 69.00000           M
1158 53.00000           M
1987 58.32055           F
712  32.79726           M
```

## Fit gaussian mixture model for this CNVR using OSCAA
Given a CNVR, simultaneously detect CNVs for each sample and  estimate the CNV-disease association adjusting age and sex. Assume only deletions of CNVs are associated with disease risk.
```
CNVR.res<-Known.region.test(position   = example.data$position,
                            signal.mat = example.data$example.lrr,
                            phenotype  = example.data$clinic.info$disease,
                            phenotype.type  = "binary",
                            dt.nCNV    = 1,
                            label.initial = "GMM",
                            Dim.reduction = "PCA",
                            assumption    = "deletion",
                            max.clusters  = 5,
                            covariates    = example.data$clinic.info[,c("Gender","age")],
                            covariates.signal.names = c("Gender"),
                            covariates.pheno.names  = c("Gender","age"),
                            overall.alpha=0.05,
                            applyldf = FALSE,
                            smooth   = FALSE,
                            max.iter = 100,
                            tol      = 0.01)
```
## Cluster plot
PCA plot with each sample labeled with identified CNV
```
OSCAA.plot(OSA.PCA.res=CNVR.res[["sig.CNR.test"]][[1]])
```

![cluster plot](https://github.com/FeifeiXiao-lab/OSCAA/blob/ee2bf75802f9dc013781ec3ef212ba040dc8e30a/image/cluster.plot.png)

## Example of phenotype models
Three hypotheses are provided in the phenotype model: Model 1 is fitted for the hypothesis of whether samples with abnormal copy number states are more/less likely to develop a disease.

$$ Model 1: E(logit(y=1))= \beta_{0} + beta_{1} * is.CNV $$

where y denotes the binary disease status while is.CNV = 1 for samples with a CNV and 0 otherwise. Model 2 is fitted for the hypothesis of whether sample with deletions of copy numbers are more/less likely to develop a disease.

$$ Model 2: E(logit(y=1))= \beta_{0} + beta_{1} * is.deletion $$

where is.deletion = 1 for samples with a deletions of copy numbers and 0 otherwise. Model 3 is fitted for the hypothesis of linear relationship between disease risk and copy number.

$$ Model 3: E(logit(y=1))= \beta_{0} + beta_{1} * CN $$

where CN represents the copy numbers.


