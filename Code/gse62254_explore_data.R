#' ---
#' title: Gastric Cancer Datasets - GSE62254
#' author: Mik Black
#' date: 26 March 2018
#' output: github_document
#' ---

#' Load required packages
library(here)
library(dplyr)
library(ggplot2)

#' Set output wider (makes things look a bit nicer: topTable etc)
options(width=80)

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

#' The sample names should be included with the clinical data, but the samples are 
#' in the same order as in the expresion data, so we can add those as a new variable:
gse62254_clinDat = gse62254_clinDat %>% mutate(., PatID = colnames(gse62254_expDat))

#' Now we have patients IDs in a column called PatID: 
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

#' Extracting expression information for a single gene (e.g., Androgen Receptor, AR).
ar_data = gse62254_expDat[rownames(gse62254_expDat) == "AR"]

#' Basic histogram of the AR data
ar_data %>% as.data.frame() %>%  
  ggplot(aes(x=.)) + 
  geom_histogram() + 
  ggtitle("Androgran receptor (AR) expression") +
  xlab("Log2 expression")

#' Boxplot of AR expression versus Molecular Subtype
cbind(ar_data, as.factor(molSub)) %>% as.data.frame() %>%  
  ggplot(aes(x=molSub, y=ar_data, group=molSub, colour=molSub)) + 
  geom_boxplot() + 
  ggtitle("Androgen receptor (AR) expression versus molecular subtype") +
  xlab("Molecular Subtype") +
  ylab("Log2 expression")

#' Differential expression analysis according to Lauren Classification.
#' <BR>

#' Lauren Classification has three levels:
table(lc)

#' Let's get rid of the "Mixed" samples
mixed = which(lc=="Mixed")

#' These are the indexes for the mixed samples
mixed

#' using the minus sign we can remove these from the Lauren data:
lc[-mixed] %>%  table()

#' Since 
#' `lc` is a factor, it still knows about the Mixed class, but it correctly 
#' reports that there are none present.  If we want to remove that class entirely,
#' we can convert lc to a vector:
lc[-mixed] %>%  as.vector() %>%  table()

#' Use this for our analysis:
lc_no_mixed = lc[-mixed] %>%  as.vector()
gse62254_expDat_no_mixed = gse62254_expDat[,-mixed]

#' Set up for limma analysis
library(limma)

#' Create design matrix with Mixed samples excluded
design = model.matrix(~lc_no_mixed)
#' Check first few rows
head(design)
#' Second column should match table created above
table(design[,2])

#' Should now be able to use `lc_no_mixed` and `deisgn` to run `limma` analysis
#' and detect genes that are differentially expressed between the Diffuse and 
#' Intestinal classes.
