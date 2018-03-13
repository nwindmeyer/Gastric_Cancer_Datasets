#' ---
#' title: Gastric Cancer Datasets - GSE62254
#' author: Mik Black
#' date: 13 March 2018
#' output: github_document
#' ---

library(here)

#' Load GSE62254 Gastric Cancer data set
load(here('Data/gse62254_gastric_cancer.RData'))

#' What was loaded?
ls()

#' What sort of objects are these?
class(gse62254_clinDat)
class(gse62254_expDat)

#' Dimensions
dim(gse62254_clinDat)
dim(gse62254_expDat)

#' Variable names for clinical data set
names(gse62254_clinDat)

#' - dfsEvent: Disease Free Survival Event (0,1)
#' - dfsTime: Time of Disease Free Survival Event (months)
#' - lc: Lauren Classification
#' - molSub: Molecular Subtype (Cristescu et al)
#' - stage: tumour stage
#' - region: anatomic region of tumour
#' - gender: patient gender
#' - ageCat: patient age category


#' First few lines of the clinical data:
head(gse62254_clinDat)

#' First 5 rows and columns of the gene expression data: rows are genes (row names 
#' are gene symbols) and columns are tumours
gse62254_expDat[1:5,1:5]

#' Table of Lauren Class data
table(gse62254_clinDat$lc)

#' Can "attach" clinical data to make the variables more accessible:
attach(gse62254_clinDat)

#' Now the variables can be referenced directly
table(lc)

#' Table of Lauren Class versus Molecular Subtype
table(lc, molSub)


