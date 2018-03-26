#' ---
#' title: Gastric Cancer - Tom's limma code
#' author: Mik Black
#' date: 27 March 2018
#' output: github_document
#' ---

#' Mik's setup:

library(here)
library(limma)

load(here('Data/gse62254_gastric_cancer.RData'))

clinData = gse62254_clinDat
attach(clinData)

#' Tom's code:

##every intestinal tumour sample ####
intest_rownames <- rownames(clinData[clinData$lc == 'Intestinal',])
##creates an empty dataframe
intest_data_frame <- data.frame(c(1:20515))
##for every column in expdata
for(i in 1:ncol(gse62254_expDat)){
  ##if that column is intestinal sample
  if(i %in% intest_rownames == TRUE ){
    ##add it to the dataframe
    intest_data_frame <- cbind(intest_data_frame, gse62254_expDat[,i])
  }
}

##every intestinal tumour sample ####
diffuse_rownames <- rownames(clinData[clinData$lc == 'Diffuse',])
##creates an empty dataframe
diffuse_data_frame <- data.frame(c(1:20515))
##for every column in expdata
for(i in 1:ncol(gse62254_expDat)){
  ##if that column is intestinal sample
  if(i %in% diffuse_rownames == TRUE ){
    ##add it to the dataframe
    diffuse_data_frame <- cbind(diffuse_data_frame, gse62254_expDat[,i])
  }
}

#' First issue: extra column in each data frame
table(lc)

#' Should be 142 samples
dim(diffuse_data_frame)

#' Should be 150 samples
dim(intest_data_frame)

#' Extra column at start of data frames:
head(intest_data_frame[,1:3])
head(diffuse_data_frame[,1:3])

#' Problem: This code adds a column (creates data.frame with single column of 
#' numbers from 1 to 20515):
## diffuse_data_frame <- data.frame(c(1:20515))

#' Could have fixed this here by removing the columns from each data frame, but 
#' I did it below on the combined object.

## cbind intestinal and diffuse data frames 
added_data_frame <- cbind(diffuse_data_frame, intest_data_frame)
head(added_data_frame[,1:3])

## check dimensions 
dim(added_data_frame) # 20515 genes , 293 samples (tumours) - mixed diffuse and intestinal samples removed
dim(gse62254_expDat) # 20515 genes, 300 samples (tumours) - confirmation 

#' BUT: 
dim(added_data_frame)  ##   20515 x 292 (2 extra columns)

#' 292 Diffuse (142) and Intestinal (150) samples
table(lc)


#' Remove extra columns
added_data_frame = added_data_frame[,-which(colnames(added_data_frame) == "c.1.20515.")]
dim(added_data_frame)

############# I now have a data frame that has just diffuse and intestinal GC samples 
############# matrix design 

#' NB - you had these the wrong way around: you used cbind to join diffuse (first) to 
#'      intestinal (second), but your labels were in the wrong order.<BR>
## group = rep(c("Intestinal","Diffuse"), c(150,142))

#' Corrected (switched order of labels):
group = rep(c("Diffuse","Intestinal"), c(142,150))

design = model.matrix(~ group)

## NB - the design matrix columns do not represent the two tumour types
##    - column one is the average expression levels across all tumours
##    - column two is the difference between Intestinal and Diffuse
##    - I've renamed them:
colnames(design) = c("Average", "Int_vs_Dif")

#' Now the code should work:

## Fit linear model
fit = lmFit(added_data_frame, design)
fit = eBayes(fit)
tt = topTable(fit, coef="Int_vs_Dif", adjust="BH",n=nrow(added_data_frame))
options(digits=4)
tt[1:5,]
## number of DE genes 
sum(tt$adj.P.Val<0.05)


