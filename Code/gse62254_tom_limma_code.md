Gastric Cancer - Tom's limma code
================
Mik Black
27 March 2018

Mik's setup:

``` r
library(here)
```

    ## here() starts at /Users/black/Research/Students/TomMarsland/GitHub/Gastric_Cancer_Datasets

``` r
library(limma)

load(here('Data/gse62254_gastric_cancer.RData'))

clinData = gse62254_clinDat
attach(clinData)
```

Tom's code:

``` r
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
```

First issue: extra column in each data frame

``` r
table(lc)
```

    ## lc
    ##    Diffuse Intestinal      Mixed 
    ##        142        150          8

Should be 142 samples

``` r
dim(diffuse_data_frame)
```

    ## [1] 20515   143

Should be 150 samples

``` r
dim(intest_data_frame)
```

    ## [1] 20515   151

Extra column at start of data frames:

``` r
head(intest_data_frame[,1:3])
```

    ##          c.1.20515. gse62254_expDat[, i] gse62254_expDat[, i].1
    ## A1BG              1             1.628701              1.6141138
    ## A1BG-AS1          2             1.494905              1.4030464
    ## A1CF              3             1.671671              2.3603470
    ## A2M               4             2.696250              3.1825163
    ## A2M-AS1           5             1.153275              1.1918744
    ## A2ML1             6             1.089608              0.9922849

``` r
head(diffuse_data_frame[,1:3])
```

    ##          c.1.20515. gse62254_expDat[, i] gse62254_expDat[, i].1
    ## A1BG              1            1.6495550               1.671960
    ## A1BG-AS1          2            1.4075329               1.514318
    ## A1CF              3            2.5890791               1.517690
    ## A2M               4            3.2541463               2.844064
    ## A2M-AS1           5            1.2974992               1.125547
    ## A2ML1             6            0.9948127               1.064056

Problem: This code adds a column (creates data.frame with single column of numbers from 1 to 20515):

``` r
## diffuse_data_frame <- data.frame(c(1:20515))
```

Could have fixed this here by removing the columns from each data frame, but I did it below on the combined object.

``` r
## cbind intestinal and diffuse data frames 
added_data_frame <- cbind(diffuse_data_frame, intest_data_frame)
head(added_data_frame[,1:3])
```

    ##          c.1.20515. gse62254_expDat[, i] gse62254_expDat[, i].1
    ## A1BG              1            1.6495550               1.671960
    ## A1BG-AS1          2            1.4075329               1.514318
    ## A1CF              3            2.5890791               1.517690
    ## A2M               4            3.2541463               2.844064
    ## A2M-AS1           5            1.2974992               1.125547
    ## A2ML1             6            0.9948127               1.064056

``` r
## check dimensions 
dim(added_data_frame) # 20515 genes , 293 samples (tumours) - mixed diffuse and intestinal samples removed
```

    ## [1] 20515   294

``` r
dim(gse62254_expDat) # 20515 genes, 300 samples (tumours) - confirmation 
```

    ## [1] 20515   300

BUT:

``` r
dim(added_data_frame)  ##   20515 x 292 (2 extra columns)
```

    ## [1] 20515   294

292 Diffuse (142) and Intestinal (150) samples

``` r
table(lc)
```

    ## lc
    ##    Diffuse Intestinal      Mixed 
    ##        142        150          8

Remove extra columns

``` r
added_data_frame = added_data_frame[,-which(colnames(added_data_frame) == "c.1.20515.")]
dim(added_data_frame)
```

    ## [1] 20515   292

``` r
############# I now have a data frame that has just diffuse and intestinal GC samples 
############# matrix design 
```

NB - you had these the wrong way around: you used cbind to join diffuse (first) to intestinal (second), but your labels were in the wrong order.<BR>

``` r
## group = rep(c("Intestinal","Diffuse"), c(150,142))
```

Corrected (switched order of labels):

``` r
group = rep(c("Diffuse","Intestinal"), c(142,150))

design = model.matrix(~ group)

## NB - the design matrix columns do not represent the two tumour types
##    - column one is the average expression levels across all tumours
##    - column two is the difference between Intestinal and Diffuse
##    - I've renamed them:
colnames(design) = c("Average", "Int_vs_Dif")
```

Now the code should work:

``` r
## Fit linear model
fit = lmFit(added_data_frame, design)
fit = eBayes(fit)
tt = topTable(fit, coef="Int_vs_Dif", adjust="BH",n=nrow(added_data_frame))
options(digits=4)
tt[1:5,]
```

    ##             logFC AveExpr      t   P.Value adj.P.Val     B
    ## IL6ST     -0.1640   3.329 -8.951 4.174e-17 7.186e-13 28.24
    ## MS4A7     -0.2836   2.580 -8.852 8.422e-17 7.186e-13 27.55
    ## GYPC      -0.1908   2.242 -8.821 1.051e-16 7.186e-13 27.34
    ## LY86      -0.2166   2.197 -8.617 4.390e-16 2.252e-12 25.94
    ## C10orf128 -0.2283   2.095 -8.419 1.728e-15 7.090e-12 24.60

``` r
## number of DE genes 
sum(tt$adj.P.Val<0.05)
```

    ## [1] 7115
