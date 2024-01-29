GACT and QGG Installation, Database Download, and Overview
================
Peter Sørensen
2024-01-29

``` r
# Install new version of gact and qgg
# Make sure you install in a new R session
# library(devtools)
# devtools::install_github("psoerensen/gact")
# devtools::install_github("psoerensen/qgg")


# Download and install gact database
library(gact)
dbdir <- "C:/Users/au223366/Dropbox/Projects/balder/gdtdb"
#GAlist <- gact(version="hsa.0.0.1", dbdir=dbdir, task="download")
#saveRDS(GAlist, file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/GAlist_hsa.0.0.1.rds")


# Overview of database content

GAlist <- readRDS(file="C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/GAlist_hsa.0.0.1.rds")

## GWAS Summary Statistics
GAlist$studies
```

    ##          id      file trait   type  gender      n ncase ncontrol     neff
    ## GWAS1 GWAS1 GWAS1.txt  T2DM binary    both 456236 55927   400309 49071.27
    ## GWAS2 GWAS2 GWAS2.txt   CAD binary    both 184305 60801   123504 40743.15
    ## GWAS3 GWAS3 GWAS3.txt  T2DM binary    both 933970 80154   853816 73275.12
    ## GWAS4 GWAS4 GWAS4.txt  T2DM binary   males 425662 41846   383816 37732.20
    ## GWAS5 GWAS5 GWAS5.txt  T2DM binary females 464379 30053   434326 28108.07
    ## GWAS6 GWAS6 GWAS6.txt  T2DM binary    both 898130 74124   824006 68006.44
    ## GWAS7 GWAS7 GWAS7.txt  T2DM binary    both 283423 56268   227155 45097.11
    ## GWAS8 GWAS8 GWAS8.txt  T2DM binary    both  49492 16540    32952 11012.41
    ##           reference                                        source
    ## GWAS1 PMID:30297969 Mahajan.NatGenet2018b.T2D-noUKBB.European.txt
    ## GWAS2 PMID:26343387                      CARDIoGRAMplusC4D.txt.gz
    ## GWAS3 PMID:35551307                   DIAMANTE-EUR.sumstat.txt.gz
    ## GWAS4 PMID:30297969   Mahajan.NatGenet2018b.T2D.MALE.European.txt
    ## GWAS5 PMID:30297969 Mahajan.NatGenet2018b.T2D.FEMALE.European.txt
    ## GWAS6 PMID:30297969        Mahajan.NatGenet2018b.T2D.European.txt
    ## GWAS7 PMID:35551307                   DIAMANTE-EAS.sumstat.txt.gz
    ## GWAS8 PMID:35551307                   DIAMANTE-SAS.sumstat.txt.gz

``` r
## Accessing the 'dirs' Slot
GAlist$dirs
```

    ##                                                                glist 
    ##    "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/glist" 
    ##                                                                gstat 
    ##    "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gstat" 
    ##                                                                gsets 
    ##    "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gsets" 
    ##                                                                 gsea 
    ##     "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gsea" 
    ##                                                               gbayes 
    ##   "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gbayes" 
    ##                                                                 gtex 
    ##     "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gtex" 
    ##                                                                 gwas 
    ##     "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gwas" 
    ##                                                                 ldsc 
    ##     "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/ldsc" 
    ##                                                               marker 
    ##   "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/marker" 
    ##                                                               drugdb 
    ##   "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/drugdb" 
    ##                                                             download 
    ## "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/download"

``` r
## Listing Directories and Files
### All Directories Recursively
list.dirs(GAlist$dirs, recursive = TRUE)
```

    ##  [1] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/download"                     
    ##  [2] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/drugdb"                       
    ##  [3] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gbayes"                       
    ##  [4] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/glist"                        
    ##  [5] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gsea"                         
    ##  [6] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gsets"                        
    ##  [7] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gstat"                        
    ##  [8] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gtex"                         
    ##  [9] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gtex/GTEx_Analysis_v7_eQTL"   
    ## [10] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gtex/GTEx_Analysis_v8_eQTL"   
    ## [11] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/gwas"                         
    ## [12] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/ldsc"                         
    ## [13] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/marker"                       
    ## [14] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/marker/1000G_EUR_Phase3_plink"
    ## [15] "C:/Users/au223366/Dropbox/Projects/balder/gdtdb/hsa.0.0.1/marker/ukb"

``` r
### Files under 'marker' Directory
list.files(GAlist$dirs["marker"])
```

    ##  [1] "1000G_EUR_Phase3_plink"            "1000G_EUR_Phase3_plink.zip"       
    ##  [3] "g1000_afr.bed"                     "g1000_afr.bim"                    
    ##  [5] "g1000_afr.fam"                     "g1000_afr.synonyms"               
    ##  [7] "g1000_afr.zip"                     "g1000_amr.bed"                    
    ##  [9] "g1000_amr.bim"                     "g1000_amr.fam"                    
    ## [11] "g1000_amr.synonyms"                "g1000_amr.zip"                    
    ## [13] "g1000_eas.bed"                     "g1000_eas.bim"                    
    ## [15] "g1000_eas.fam"                     "g1000_eas.synonyms"               
    ## [17] "g1000_eas.zip"                     "g1000_eas_filtered.bed"           
    ## [19] "g1000_eas_filtered.bim"            "g1000_eas_filtered.fam"           
    ## [21] "g1000_eur.bed"                     "g1000_eur.bim"                    
    ## [23] "g1000_eur.fam"                     "g1000_eur.synonyms"               
    ## [25] "g1000_eur.zip"                     "g1000_eur_filtered.bed"           
    ## [27] "g1000_eur_filtered.bim"            "g1000_eur_filtered.fam"           
    ## [29] "g1000_sas.bed"                     "g1000_sas.bim"                    
    ## [31] "g1000_sas.fam"                     "g1000_sas.synonyms"               
    ## [33] "g1000_sas.zip"                     "g1000_sas_filtered.bed"           
    ## [35] "g1000_sas_filtered.bim"            "g1000_sas_filtered.fam"           
    ## [37] "Glist_1000G.rds"                   "Glist_1000G_eas_filtered.rds"     
    ## [39] "Glist_1000G_eur_filtered.rds"      "Glist_1000G_sas_filtered.rds"     
    ## [41] "markers.txt.gz"                    "markers_1000G.txt.gz"             
    ## [43] "markers_1000G_eas_filtered.txt.gz" "markers_1000G_eur_filtered.txt.gz"
    ## [45] "markers_1000G_sas_filtered.txt.gz" "README"                           
    ## [47] "ukb"

``` r
### Files under 'gstat' Directory
list.files(GAlist$dirs["gstat"])
```

    ## [1] "GWAS_information.csv" "GWAS1.txt.gz"         "GWAS2.txt.gz"        
    ## [4] "GWAS3.txt.gz"         "GWAS4.txt.gz"         "GWAS5.txt.gz"        
    ## [7] "GWAS6.txt.gz"         "GWAS7.txt.gz"         "GWAS8.txt.gz"

``` r
### Files under 'gsea' Directory
list.files(GAlist$dirs["gsea"])
```

    ##  [1] "GWAS1_vegas.rds"               "GWAS10_vegas.rds"             
    ##  [3] "GWAS11_vegas.rds"              "GWAS12_vegas.rds"             
    ##  [5] "GWAS13_vegas.rds"              "GWAS14_vegas.rds"             
    ##  [7] "GWAS15_vegas.rds"              "GWAS16_vegas.rds"             
    ##  [9] "GWAS17_vegas.rds"              "GWAS18_vegas.rds"             
    ## [11] "GWAS2_vegas.rds"               "GWAS22_vegas.rds"             
    ## [13] "GWAS25_vegas.rds"              "GWAS3_vegas.rds"              
    ## [15] "GWAS4_vegas.rds"               "GWAS5_vegas.rds"              
    ## [17] "GWAS6_vegas.rds"               "GWAS7_vegas.rds"              
    ## [19] "GWAS8_vegas.rds"               "p_vegas.csv"                  
    ## [21] "P_vegas.rds"                   "p_vegas.txt.gz"               
    ## [23] "P_vegas_filtered.rds"          "res_vegas_eas_filtered.rds"   
    ## [25] "res_vegas_eur_filtered.rds"    "res_vegas_sas_filtered.rds"   
    ## [27] "summary_vegas_filtered.csv"    "summary_vegas_filtered.txt.gz"
    ## [29] "Z_vegas.rds"                   "z_vegas.txt.gz"               
    ## [31] "Z_vegas_filtered.rds"
