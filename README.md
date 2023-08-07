# OSCAA
A Two-Dimensional Gaussian Mixture Model for Copy Number Variation Association Analysis

## Author
Xuanxuan Yu, Xizhi Luo, Guoshuai Cai, Feifei Xiao

## Description
The One-Stage CNV-disease Association Analysis (OSCAA) method utilizes a two-dimensional GMM model to assess the association between CNV and disease risk in a known CNVR with a one-stage framework. The two-dimensional GMM model is factorized into a signal model, a phenotype model, and a copy number model. To better differentiate samples with different copy numbers, we model the joint distribution of the first two PCs for samples with the same copy number in the signal model. The phenotype model comprises a generalized linear model (GLM) or a logistic model fitted for continuous or binary phenotype measurements, respectively. In addition, the phenotype model offers two hypotheses: the hypothesis of linear relationship between disease risk and copy number, and the hypothesis of whether samples with abnormal copy number states are more likely to develop a disease. The copy number model estimates the proportion of each copy number state and enables the evaluation of the potential impact of covariates on the proportions. Although we demonstrate the model based on SNP array data, our method can be extended to sequencing data with appropriate data normalization procedures.

## install OSCAA
```
library(devtools)
install_github("....")
```

# data
```
library("OSCAA")
data(exmaple.data)
dim(example.data$example.lrr)
[1] 2829   15
example.data$position
```
