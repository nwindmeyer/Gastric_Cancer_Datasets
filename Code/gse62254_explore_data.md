Gastric Cancer Datasets - GSE62254
================
Mik Black
13 March 2018

``` r
library(here)
```

    ## here() starts at /Users/black/Research/Students/TomMarsland/GitHub/Gastric_Cancer_Datasets

Load GSE62254 Gastric Cancer data set

``` r
load(here('Data/gse62254_gastric_cancer.RData'))
```

What was loaded?

``` r
ls()
```

    ## [1] "gse62254_clinDat" "gse62254_expDat"

What sort of objects are these?

``` r
class(gse62254_clinDat)
```

    ## [1] "data.frame"

``` r
class(gse62254_expDat)
```

    ## [1] "matrix"

Dimensions

``` r
dim(gse62254_clinDat)
```

    ## [1] 300   8

``` r
dim(gse62254_expDat)
```

    ## [1] 20515   300

Variable names for clinical data set

``` r
names(gse62254_clinDat)
```

    ## [1] "dfsEvent" "dfsTime"  "lc"       "molSub"   "stage"    "region"  
    ## [7] "gender"   "ageCat"

-   dfsEvent: Disease Free Survival Event (0,1)
-   dfsTime: Time of Disease Free Survival Event (months)
-   lc: Lauren Classification
-   molSub: Molecular Subtype (Cristescu et al)
-   stage: tumour stage
-   region: anatomic region of tumour
-   gender: patient gender
-   ageCat: patient age category First few lines of the clinical data:

``` r
head(gse62254_clinDat)
```

    ##   dfsEvent dfsTime         lc    molSub stage region gender ageCat
    ## 1       NA    3.97 Intestinal       MSI     2   body   Male  65-69
    ## 2       NA    4.03 Intestinal       MSI     2   body Female  65-69
    ## 3        0   74.97    Diffuse MSS/TP53+     2 antrum Female    <55
    ## 4        0   89.77    Diffuse       MSI     2 antrum   Male  65-69
    ## 5        0   84.60    Diffuse MSS/TP53-     3 antrum   Male  65-69
    ## 6        1    5.77      Mixed MSS/TP53-     2 antrum   Male  55-64

First 5 rows and columns of the gene expression data: rows are genes (row names are gene symbols) and columns are tumours

``` r
gse62254_expDat[1:5,1:5]
```

    ##          GSM1523727 GSM1523728 GSM1523729 GSM1523744 GSM1523745
    ## A1BG       1.628701   1.614114   1.649555   1.671960   1.636245
    ## A1BG-AS1   1.494905   1.403046   1.407533   1.514318   1.471773
    ## A1CF       1.671671   2.360347   2.589079   1.517690   2.335692
    ## A2M        2.696250   3.182516   3.254146   2.844064   3.258593
    ## A2M-AS1    1.153275   1.191874   1.297499   1.125547   1.145582

Table of Lauren Class data

``` r
table(gse62254_clinDat$lc)
```

    ## 
    ##    Diffuse Intestinal      Mixed 
    ##        142        150          8

Can "attach" clinical data to make the variables more accessible:

``` r
attach(gse62254_clinDat)
```

Now the variables can be referenced directly

``` r
table(lc)
```

    ## lc
    ##    Diffuse Intestinal      Mixed 
    ##        142        150          8

Table of Lauren Class versus Molecular Subtype

``` r
table(lc, molSub)
```

    ##             molSub
    ## lc           EMT MSI MSS/TP53- MSS/TP53+
    ##   Diffuse     38  20        46        38
    ##   Intestinal   8  43        59        40
    ##   Mixed        0   5         2         1
