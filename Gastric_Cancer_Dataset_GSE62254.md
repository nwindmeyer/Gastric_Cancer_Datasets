Gastric Cancer Dataset Practise - GSE62254
================
Natasha Windmeyer
19/03/2021

This Practise document is using the Affymetrix microarray dataset
(GSE62254) - re-normalised using RMA without background
    correction

``` r
install.packages('here')
```

``` r
library(here)
```

    ## here() starts at /Users/natashawindmeyer/Documents/University MSc/R Data/Gastric_Cancer_Datasets

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
```

The ‘here()’ starts at /Users/black/GitHub/Gastric\_Cancer\_Datasets

``` r
options(width=80)
```

This sets the output wider to 80, making things look a bit nicer

``` r
load(here('Data/gse62254_gastric_cancer.RData'))
```

This loads the GSE62254 Gastric Cancer Dataset

``` r
ls()
```

    ## [1] "gse62254_clinDat" "gse62254_expDat"

This confirms what was loaded

``` r
class(gse62254_clinDat)
```

    ## [1] "data.frame"

``` r
class(gse62254_expDat)
```

    ## [1] "matrix" "array"

This tells you what type of objects the datasets are

``` r
dim(gse62254_clinDat)
```

    ## [1] 300   8

``` r
dim(gse62254_expDat)
```

    ## [1] 20184   300

This gives you the dimensions of the datasets

The Variable names for the Clinical
    Dataset

``` r
names(gse62254_clinDat)
```

    ## [1] "dfsEvent" "dfsTime"  "lc"       "molSub"   "stage"    "region"   "gender"  
    ## [8] "ageCat"

dfsEvent: Disease Free Survival Event (0,1) dfsTime: Time of Disease
Free Survival Event (months) lc: Lauren Classification molSub: Molecular
Subtype (Cristescu et al) stage: Tumour stage region: Anatomic region of
tumour gender: Patient gender ageCat: Patient age
    category

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

Gives you the first 6 lines of the clinical dataset

The sample names are not included in the clinical data, so we need to
add this as a new
column

``` r
gse62254_clinDat = gse62254_clinDat %>% mutate(, PatID = colnames(gse62254_expDat))
head(gse62254_clinDat)
```

    ##   dfsEvent dfsTime         lc    molSub stage region gender ageCat      PatID
    ## 1       NA    3.97 Intestinal       MSI     2   body   Male  65-69 GSM1523727
    ## 2       NA    4.03 Intestinal       MSI     2   body Female  65-69 GSM1523728
    ## 3        0   74.97    Diffuse MSS/TP53+     2 antrum Female    <55 GSM1523729
    ## 4        0   89.77    Diffuse       MSI     2 antrum   Male  65-69 GSM1523744
    ## 5        0   84.60    Diffuse MSS/TP53-     3 antrum   Male  65-69 GSM1523745
    ## 6        1    5.77      Mixed MSS/TP53-     2 antrum   Male  55-64 GSM1523746
