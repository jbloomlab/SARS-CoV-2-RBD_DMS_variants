Collapse barcodes to final per-RBD/mutant phenotype scores
================
Tyler Starr
07/28/2021

-   [Setup](#setup)
-   [Calculate per-variant mean scores within
    replicates](#calculate-per-variant-mean-scores-within-replicates)
-   [Calculate per-mutant score across
    libraries](#calculate-per-mutant-score-across-libraries)
-   [Correlations among backgrounds and to prior Wuhan-Hu-1 DMS
    data](#correlations-among-backgrounds-and-to-prior-wuhan-hu-1-dms-data)
-   [Heatmaps!](#heatmaps)

This notebook reads in the per-barcode titration Kds and expression
measurements from the `compute_binding_Kd` and
`compute_expression_meanF` scripts. It synthesizes these two sets of
results and calculates the final ‘mean’ phenotypes for each variant, and
generates some coverage and QC analyses.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","egg")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
invisible(lapply(packages, library, character.only=T))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$final_variant_scores_dir)){
  dir.create(file.path(config$final_variant_scores_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.7-GCC-8.3.0/lib/libopenblas_haswellp-r0.3.7.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] egg_0.4.5         gridExtra_2.3     forcats_0.4.0     stringr_1.4.0    
    ##  [5] dplyr_0.8.3       purrr_0.3.3       readr_1.3.1       tidyr_1.0.0      
    ##  [9] tibble_3.0.2      ggplot2_3.3.0     tidyverse_1.3.0   data.table_1.12.8
    ## [13] yaml_2.2.0        knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0 xfun_0.11        haven_2.2.0      colorspace_1.4-1
    ##  [5] vctrs_0.3.1      generics_0.0.2   htmltools_0.4.0  rlang_0.4.7     
    ##  [9] pillar_1.4.5     glue_1.3.1       withr_2.1.2      DBI_1.1.0       
    ## [13] dbplyr_1.4.2     modelr_0.1.5     readxl_1.3.1     lifecycle_0.2.0 
    ## [17] munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0 rvest_0.3.5     
    ## [21] evaluate_0.14    fansi_0.4.0      broom_0.7.0      Rcpp_1.0.3      
    ## [25] scales_1.1.0     backports_1.1.5  jsonlite_1.6     fs_1.3.1        
    ## [29] hms_0.5.2        digest_0.6.23    stringi_1.4.3    grid_3.6.2      
    ## [33] cli_2.0.0        tools_3.6.2      magrittr_1.5     crayon_1.3.4    
    ## [37] pkgconfig_2.0.3  ellipsis_0.3.0   xml2_1.2.2       reprex_0.3.0    
    ## [41] lubridate_1.7.4  assertthat_0.2.1 rmarkdown_2.0    httr_1.4.1      
    ## [45] rstudioapi_0.10  R6_2.4.1         compiler_3.6.2

## Setup

Read in tables of per-barcode expression and binding Kd measurements and
combine.

``` r
dt_bind <- data.table(read.csv(config$Titeseq_Kds_file),stringsAsFactors=F)
dt_expr <- data.table(read.csv(config$expression_sortseq_file),stringsAsFactors=F)
```

## Calculate per-variant mean scores within replicates

Calculate the mean binding and expression score collapsed by genotype.
Also output the number of barcodes across which a variant score was
determined in each library.

``` r
dt_bind[is.na(log10Ka),TiteSeq_avgcount:=NA]
dt_expr[is.na(expression),expr_count:=NA]

dt_bind[,mean_bind:=mean(log10Ka,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_bind[,sd_bind:=sd(log10Ka,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_bind[,n_bc_bind:=sum(!is.na(log10Ka)),by=c("library","target","variant_class","aa_substitutions")]
dt_bind[,avg_count_bind:=mean(TiteSeq_avgcount,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]

dt_bind <- unique(dt_bind[,.(library,target,variant_class,aa_substitutions,n_aa_substitutions,mean_bind,sd_bind,n_bc_bind,avg_count_bind)])

dt_expr[,mean_expr:=mean(expression,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_expr[,sd_expr:=sd(expression,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]
dt_expr[,n_bc_expr:=sum(!is.na(expression)),by=c("library","target","variant_class","aa_substitutions")]
dt_expr[,avg_count_expr:=mean(expr_count,na.rm=T),by=c("library","target","variant_class","aa_substitutions")]

dt_expr <- unique(dt_expr[,.(library,target,variant_class,aa_substitutions,n_aa_substitutions,mean_expr,sd_expr,n_bc_expr,avg_count_expr)])
```

Some QC plots. First, look at distribution of number barcodes for
binding and expression measurements for single mutant detemrinations.
These are ‘left-justified’ histograms, so the leftmost bar represents
the number of genotypes for which no barcodes were collapsed to final
measurement in a pool.

``` r
par(mfrow=c(2,3))
hist(dt_bind[library=="pool1A" & variant_class=="1 nonsynonymous",n_bc_bind],main="pool1A, bind",right=F,breaks=max(dt_bind[library=="pool1A" & variant_class=="1 nonsynonymous",n_bc_bind],na.rm=T),xlab="")
hist(dt_bind[library=="pool1B" & variant_class=="1 nonsynonymous",n_bc_bind],main="pool1B, bind",right=F,breaks=max(dt_bind[library=="pool1B" & variant_class=="1 nonsynonymous",n_bc_bind],na.rm=T),xlab="")
hist(dt_bind[library=="pool2A" & variant_class=="1 nonsynonymous",n_bc_bind],main="pool2A, bind",right=F,breaks=max(dt_bind[library=="pool2A" & variant_class=="1 nonsynonymous",n_bc_bind],na.rm=T),xlab="")
hist(dt_expr[library=="pool1" & variant_class=="1 nonsynonymous",n_bc_expr],main="pool1, expr",right=F,breaks=max(dt_expr[library=="pool1" & variant_class=="1 nonsynonymous",n_bc_expr],na.rm=T),xlab="number barcodes collapsed")
plot(0,type='n',axes=FALSE,ann=F)
hist(dt_expr[library=="pool2" & variant_class=="1 nonsynonymous",n_bc_expr],main="pool2, expr",right=F,breaks=max(dt_expr[library=="pool2" & variant_class=="1 nonsynonymous",n_bc_expr],na.rm=T),xlab="number barcodes collapsed")
```

<img src="collapse_scores_files/figure-gfm/hist_n_bc_per_mutant-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/histogram_n_bc_per_geno_sep-libs.pdf",sep=""),useDingbats=F))
```

What about how SEM tracks with number of barcodes collapsed? This could
help for choosing a minimum number of barcodes to use.

``` r
par(mfrow=c(2,3))
plot(dt_bind[library=="pool1A" & variant_class=="1 nonsynonymous",n_bc_bind],
     dt_bind[library=="pool1A" & variant_class=="1 nonsynonymous",sd_bind/sqrt(n_bc_bind)],
     pch=19,col="#00000005",main="pool1A, bind",ylab="SEM",xlab="number barcodes collapsed")
plot(dt_bind[library=="pool1B" & variant_class=="1 nonsynonymous",n_bc_bind],
     dt_bind[library=="pool1B" & variant_class=="1 nonsynonymous",sd_bind/sqrt(n_bc_bind)],
     pch=19,col="#00000005",main="pool1B, bind",ylab="SEM",xlab="number barcodes collapsed")
plot(dt_bind[library=="pool2A" & variant_class=="1 nonsynonymous",n_bc_bind],
     dt_bind[library=="pool2A" & variant_class=="1 nonsynonymous",sd_bind/sqrt(n_bc_bind)],
     pch=19,col="#00000005",main="pool2A, bind",ylab="SEM",xlab="number barcodes collapsed")
plot(dt_expr[library=="pool1" & variant_class=="1 nonsynonymous",n_bc_expr],
     dt_expr[library=="pool1" & variant_class=="1 nonsynonymous",sd_expr/sqrt(n_bc_expr)],
     pch=19,col="#00000005",main="pool1, expr",ylab="SEM",xlab="number barcodes collapsed")
plot(0,type='n',axes=FALSE,ann=F)
plot(dt_expr[library=="pool2" & variant_class=="1 nonsynonymous",n_bc_expr],
     dt_expr[library=="pool2" & variant_class=="1 nonsynonymous",sd_expr/sqrt(n_bc_expr)],
     pch=19,col="#00000005",main="pool2, expr",ylab="SEM",xlab="number barcodes collapsed")
```

<img src="collapse_scores_files/figure-gfm/sem_v_n-bc-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/sem_v_n-bc.pdf",sep=""),useDingbats=F))
```

Format into a ‘mutation lookup table’, where we focus just on the single
mutants (and wildtype), breakup the string of mutations, and fill in the
table to also include any missing mutants.

``` r
dt_mutant_bind <- dt_bind[variant_class %in% "1 nonsynonymous",]

#split mutation string
#define function to apply
split_mut <- function(x){
  split <- strsplit(x,split="")[[1]]
  return(list(split[1],as.numeric(paste(split[2:(length(split)-1)],collapse="")),split[length(split)]))
}
dt_mutant_bind[,c("wildtype","position","mutant"):=split_mut(as.character(aa_substitutions)),by=aa_substitutions]

dt_mutant_bind <- dt_mutant_bind[,.(library,target,wildtype,position,mutant,mean_bind,sd_bind,n_bc_bind,avg_count_bind)]

aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
#fill out missing values in table with a hideous loop, so the table is complete for all mutaitons (including those that are missing). If you are somebody who is reading this code, I apologize.
for(lib in c("pool1A","pool1B","pool2A")){
  for(bg in as.character(unique(dt_mutant_bind$target))){
    for(pos in 1:max(dt_mutant_bind$position)){
      for(aa in aas){
        if(!(aa %in% as.character(dt_mutant_bind[library==lib & target==bg & position==pos,mutant]))){
          dt_mutant_bind <- rbind(dt_mutant_bind,list(lib, bg, dt_mutant_bind[library==lib & target==bg & position==pos,wildtype][1],pos,aa),fill=T)
        }
      }
    }
  }
}
setkey(dt_mutant_bind,library,target,position,mutant)

#fill in wildtype values -- should vectorize in data table but being so stupid so just going to write for loop
for(bg in c("Wuhan_Hu_1","E484K","N501Y","B1351")){
  for(lib in c("pool1A","pool1B","pool2A")){
    dt_mutant_bind[library==lib & target==bg & wildtype==mutant, c("mean_bind","sd_bind","n_bc_bind","avg_count_bind"):=dt_bind[library==lib & target==bg & variant_class=="wildtype",.(mean_bind,sd_bind,n_bc_bind,avg_count_bind)]]
  }
}

#add delta bind measures
for(bg in c("Wuhan_Hu_1","E484K","N501Y","B1351")){
  for(lib in c("pool1A","pool1B","pool2A")){
    ref_bind <- dt_bind[library==lib & target==bg & variant_class=="wildtype",mean_bind]
    dt_mutant_bind[library==lib & target==bg,delta_bind := mean_bind - ref_bind]
  }
}

#repeat for expr
dt_mutant_expr <- dt_expr[variant_class %in% "1 nonsynonymous",]

#split mutation string
#define function to apply
split_mut <- function(x){
  split <- strsplit(x,split="")[[1]]
  return(list(split[1],as.numeric(paste(split[2:(length(split)-1)],collapse="")),split[length(split)]))
}
dt_mutant_expr[,c("wildtype","position","mutant"):=split_mut(as.character(aa_substitutions)),by=aa_substitutions]

dt_mutant_expr <- dt_mutant_expr[,.(library,target,wildtype,position,mutant,mean_expr,sd_expr,n_bc_expr,avg_count_expr)]

aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
#fill out missing values in table with a hideous loop, so the table is complete for all mutaitons (including those that are missing). If you are somebody who is reading this code, I apologize.
for(lib in c("pool1","pool2")){
  for(bg in as.character(unique(dt_mutant_expr$target))){
    for(pos in 1:max(dt_mutant_expr$position)){
      for(aa in aas){
        if(!(aa %in% as.character(dt_mutant_expr[library==lib & target==bg & position==pos,mutant]))){
          dt_mutant_expr <- rbind(dt_mutant_expr,list(lib, bg, dt_mutant_expr[library==lib & target==bg & position==pos,wildtype][1],pos,aa),fill=T)
        }
      }
    }
  }
}
setkey(dt_mutant_expr,library,target,position,mutant)

#fill in wildtype values -- should vectorize in data table but being so stupid so just going to write for loop
for(bg in c("Wuhan_Hu_1","E484K","N501Y","B1351")){
  for(lib in c("pool1","pool2")){
    dt_mutant_expr[library==lib & target==bg & wildtype==mutant, c("mean_expr","sd_expr","n_bc_expr","avg_count_expr"):=dt_expr[library==lib & target==bg & variant_class=="wildtype",.(mean_expr,sd_expr,n_bc_expr,avg_count_expr)]]
  }
}

#add delta expr measures
for(bg in c("Wuhan_Hu_1","E484K","N501Y","B1351")){
  for(lib in c("pool1","pool2")){
    ref_expr <- dt_expr[library==lib & target==bg & variant_class=="wildtype",mean_expr]
    dt_mutant_expr[library==lib & target==bg,delta_expr := mean_expr - ref_expr]
  }
}
```

We have duplicates or triplicates for each measurement. Let’s look at
correlations!

``` r
par(mfrow=c(2,3))
x <- dt_mutant_expr[library=="pool1" & wildtype!=mutant,mean_expr]; y <- dt_mutant_expr[library=="pool2" & wildtype!=mutant,mean_expr]; plot(x,y,pch=19,col="#00000020",xlab="replicate 1",ylab="replicate 2",main="expression");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
plot(0,type='n',axes=FALSE,ann=F)
plot(0,type='n',axes=FALSE,ann=F)
x <- dt_mutant_bind[library=="pool1A" & wildtype!=mutant,mean_bind]; y <- dt_mutant_bind[library=="pool2A" & wildtype!=mutant,mean_bind]; plot(x,y,pch=19,col="#00000020",xlab="replicate 1A",ylab="replicate 2",main="binding affinity");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
x <- dt_mutant_bind[library=="pool1B" & wildtype!=mutant,mean_bind]; y <- dt_mutant_bind[library=="pool2A" & wildtype!=mutant,mean_bind]; plot(x,y,pch=19,col="#00000020",xlab="replicate 1B",ylab="replicate 2",main="binding affinity");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
x <- dt_mutant_bind[library=="pool1A" & wildtype!=mutant,mean_bind]; y <- dt_mutant_bind[library=="pool1B" & wildtype!=mutant,mean_bind]; plot(x,y,pch=19,col="#00000020",xlab="replicate 1A",ylab="replicate 1B",main="binding affinity");model <- lm(y~x);abline(model,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/replicate_correlations.pdf",sep=""),useDingbats=F))
```

## Calculate per-mutant score across libraries

Collapse down to mean from both replicates, and total n barcodes between
the two/three replicates. Also record the number of the replicates the
variant was quantified within. Note, we are currently keeping a value
even if it’s determined from a single bc fit in a single pool. Later on,
we may want to require some combination of minimum number of bcs within
or between libraries for retention.

``` r
dt_final_bind <- copy(dt_mutant_bind)

dt_final_bind[ ,bind_tot:=mean(mean_bind,na.rm=T),by=c("target","position","mutant")]
dt_final_bind[ ,delta_bind_tot:=mean(delta_bind,na.rm=T),by=c("target","position","mutant")]
dt_final_bind[ ,n_bc_bind_tot:=sum(n_bc_bind,na.rm=T),by=c("target","position","mutant")]
dt_final_bind[ ,n_libs_bind_tot:=sum(!is.na(mean_bind)),by=c("target","position","mutant")]

#switch to spike indexing of postitions
dt_final_bind$position <- dt_final_bind$position + config$site_number_offset

#add single mutation string
dt_final_bind[,mutation:=paste(wildtype,position,mutant,sep=""),by=c("wildtype","position","mutant")]

dt_final_bind <- unique(dt_final_bind[,.(target,wildtype,position,mutant,mutation,bind_tot,delta_bind_tot,n_bc_bind_tot,n_libs_bind_tot)])

#repeat expr
dt_final_expr <- copy(dt_mutant_expr)

dt_final_expr[ ,expr_tot:=mean(mean_expr,na.rm=T),by=c("target","position","mutant")]
dt_final_expr[ ,delta_expr_tot:=mean(delta_expr,na.rm=T),by=c("target","position","mutant")]
dt_final_expr[ ,n_bc_expr_tot:=sum(n_bc_expr,na.rm=T),by=c("target","position","mutant")]
dt_final_expr[ ,n_libs_expr_tot:=sum(!is.na(mean_expr)),by=c("target","position","mutant")]

#switch to spike indexing of postitions
dt_final_expr$position <- dt_final_expr$position + config$site_number_offset

#add single mutation string
dt_final_expr[,mutation:=paste(wildtype,position,mutant,sep=""),by=c("wildtype","position","mutant")]

dt_final_expr <- unique(dt_final_expr[,.(target,wildtype,position,mutant,mutation,expr_tot,delta_expr_tot,n_bc_expr_tot,n_libs_expr_tot)])

#merge together
dt_final <- merge(dt_final_bind, dt_final_expr)
setkey(dt_final,target,position,mutant)

#add the rep1 and rep2 bind and expr averages
dt_final[,bind_rep1 := dt_mutant_bind[library=="pool1A", mean_bind]]
dt_final[,bind_rep2 := dt_mutant_bind[library=="pool2A", mean_bind]]
dt_final[,bind_rep3 := dt_mutant_bind[library=="pool1B", mean_bind]]
dt_final[,expr_rep1 := dt_mutant_expr[library=="pool1", mean_expr]]
dt_final[,expr_rep2 := dt_mutant_expr[library=="pool2", mean_expr]]


#rename some of the columns
setnames(dt_final,"bind_tot","bind")
setnames(dt_final,"delta_bind_tot","delta_bind")
setnames(dt_final,"n_bc_bind_tot","n_bc_bind")
setnames(dt_final,"n_libs_bind_tot","n_libs_bind")
setnames(dt_final,"expr_tot","expr")
setnames(dt_final,"delta_expr_tot","delta_expr")
setnames(dt_final,"n_bc_expr_tot","n_bc_expr")
setnames(dt_final,"n_libs_expr_tot","n_libs_expr")
```

Merge in the Delta-bg DMS measurements, obtained from Delta repo

``` r
dt_delta <- data.table(read.csv(file=config$delta_mut_bind_expr,stringsAsFactors = F))

dt_final <- rbind(dt_final, dt_delta, fill=T)
setkey(dt_final,target,position,mutant)
```

Censor any measurements that are from \<3 bc or only sampled in a single
replicate? Don’t do this for now.

``` r
# min_bc <- 2
# min_lib <- 2
# 
# dt_final[n_bc_bind < min_bc & n_libs_bind < min_lib, c("bind","delta_bind","n_bc_bind","n_libs_bind") := list(NA,NA,NA,NA)]
# dt_final[n_bc_expr < min_bc & n_libs_expr < min_lib, c("expr","delta_expr","n_bc_expr","n_libs_expr") := list(NA,NA,NA,NA)]
```

Coverage stats on n_barcodes for different measurements in the final
pooled measurements. (Exclude delta)

``` r
par(mfrow=c(1,2))
hist(dt_final[wildtype!=mutant & target!="Delta", n_bc_bind],col="gray50",main=paste("mutant bind score,\nmedian ",median(dt_final[wildtype!=mutant & target!="Delta", n_bc_bind],na.rm=T),sep=""),right=F,breaks=max(dt_final[wildtype!=mutant & target!="Delta", n_bc_bind]),xlab="number barcodes")
hist(dt_final[wildtype!=mutant & target!="Delta", n_bc_expr],col="gray50",main=paste("mutant expr score,\nmedian ",median(dt_final[wildtype!=mutant & target!="Delta", n_bc_expr],na.rm=T),sep=""),right=F,breaks=max(dt_final[wildtype!=mutant & target!="Delta", n_bc_expr]),xlab="")
```

<img src="collapse_scores_files/figure-gfm/n_barcode_plots-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/histogram_n_bc_per_geno_pooled-libs.pdf",sep="")))
```

## Correlations among backgrounds and to prior Wuhan-Hu-1 DMS data

Look at correlations in mutation effects between each background, for
bind phenotype

``` r
par(mfrow=c(4,4))

x <- dt_final[target=="B1351",bind]; y <- dt_final[target=="Delta",bind]; plot(x,y,pch=19,col="#00000020",xlab="B1351",ylab="Delta",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="N501Y",bind]; y <- dt_final[target=="Delta",bind]; plot(x,y,pch=19,col="#00000020",xlab="N501Y",ylab="Delta",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="E484K",bind]; y <- dt_final[target=="Delta",bind]; plot(x,y,pch=19,col="#00000020",xlab="E484K",ylab="Delta",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="Wuhan_Hu_1",bind]; y <- dt_final[target=="Delta",bind]; plot(x,y,pch=19,col="#00000020",xlab="Wuhan_Hu_1",ylab="Delta",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

plot(0,type='n',axes=FALSE,ann=F)

x <- dt_final[target=="N501Y",bind]; y <- dt_final[target=="B1351",bind]; plot(x,y,pch=19,col="#00000020",xlab="N501Y",ylab="B1351",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="E484K",bind]; y <- dt_final[target=="B1351",bind]; plot(x,y,pch=19,col="#00000020",xlab="E484K",ylab="B1351",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="Wuhan_Hu_1",bind]; y <- dt_final[target=="B1351",bind]; plot(x,y,pch=19,col="#00000020",xlab="Wuhan_Hu_1",ylab="B1351",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

plot(0,type='n',axes=FALSE,ann=F)

plot(0,type='n',axes=FALSE,ann=F)

x <- dt_final[target=="E484K",bind]; y <- dt_final[target=="N501Y",bind]; plot(x,y,pch=19,col="#00000020",xlab="E484K",ylab="N501Y",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="Wuhan_Hu_1",bind]; y <- dt_final[target=="N501Y",bind]; plot(x,y,pch=19,col="#00000020",xlab="Wuhan_Hu_1",ylab="N501Y",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

plot(0,type='n',axes=FALSE,ann=F)

plot(0,type='n',axes=FALSE,ann=F)

plot(0,type='n',axes=FALSE,ann=F)

x <- dt_final[target=="Wuhan_Hu_1",bind]; y <- dt_final[target=="E484K",bind]; plot(x,y,pch=19,col="#00000020",xlab="Wuhan_Hu_1",ylab="E484K",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations_by_bg_bind-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/background_correlations_bind.pdf",sep=""),useDingbats=F))
```

Look at correlations in mutation effects between each background, for
expr phenotype

``` r
par(mfrow=c(4,4))

x <- dt_final[target=="B1351",expr]; y <- dt_final[target=="Delta",expr]; plot(x,y,pch=19,col="#00000020",xlab="B1351",ylab="Delta",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="N501Y",expr]; y <- dt_final[target=="Delta",expr]; plot(x,y,pch=19,col="#00000020",xlab="N501Y",ylab="Delta",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="E484K",expr]; y <- dt_final[target=="Delta",expr]; plot(x,y,pch=19,col="#00000020",xlab="E484K",ylab="Delta",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="Wuhan_Hu_1",expr]; y <- dt_final[target=="Delta",expr]; plot(x,y,pch=19,col="#00000020",xlab="Wuhan_Hu_1",ylab="Delta",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

plot(0,type='n',axes=FALSE,ann=F)

x <- dt_final[target=="N501Y",expr]; y <- dt_final[target=="B1351",expr]; plot(x,y,pch=19,col="#00000020",xlab="N501Y",ylab="B1351",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="E484K",expr]; y <- dt_final[target=="B1351",expr]; plot(x,y,pch=19,col="#00000020",xlab="E484K",ylab="B1351",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="Wuhan_Hu_1",expr]; y <- dt_final[target=="B1351",expr]; plot(x,y,pch=19,col="#00000020",xlab="Wuhan_Hu_1",ylab="B1351",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

plot(0,type='n',axes=FALSE,ann=F)

plot(0,type='n',axes=FALSE,ann=F)

x <- dt_final[target=="E484K",expr]; y <- dt_final[target=="N501Y",expr]; plot(x,y,pch=19,col="#00000020",xlab="E484K",ylab="N501Y",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_final[target=="Wuhan_Hu_1",expr]; y <- dt_final[target=="N501Y",expr]; plot(x,y,pch=19,col="#00000020",xlab="Wuhan_Hu_1",ylab="N501Y",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

plot(0,type='n',axes=FALSE,ann=F)

plot(0,type='n',axes=FALSE,ann=F)

plot(0,type='n',axes=FALSE,ann=F)

x <- dt_final[target=="Wuhan_Hu_1",expr]; y <- dt_final[target=="E484K",expr]; plot(x,y,pch=19,col="#00000020",xlab="Wuhan_Hu_1",ylab="E484K",main="", xlim=c(5,11.5),ylim=c(5,11.5));model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations_by_bg_expr-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/background_correlations_expr.pdf",sep=""),useDingbats=F))
```

And, look at relationship between the Wuhan-Hu-1 background, monomer
ACE2 measurements versus our previously published measurements for
binding to dimeric ACE2

``` r
dt_og <- data.table(read.csv(file=config$mut_bind_expr,stringsAsFactors = F))

dt_og$bind_new <- as.numeric(NA)
dt_og$expr_new <- as.numeric(NA)

for(i in 1:nrow(dt_og)){
  if(dt_og[i,mutant]!="*"){
    dt_og[i,"bind_new"] <- dt_final[target=="Wuhan_Hu_1" & position==dt_og[i,site_SARS2] & mutant==dt_og[i,mutant],delta_bind]
    dt_og[i,"expr_new"] <- dt_final[target=="Wuhan_Hu_1" & position==dt_og[i,site_SARS2] & mutant==dt_og[i,mutant],delta_expr]
  }
}

par(mfrow=c(1,2))
x <- dt_og[,expr_avg]; y <- dt_og[,expr_new]; plot(x,y,pch=19,col="#00000020",xlab="original DMS",ylab="new DMS",main="expression");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")

x <- dt_og[,bind_avg]; y <- dt_og[,bind_new]; plot(x,y,pch=19,col="#00000020",xlab="original DMS (ACE2 dimer)",ylab="new DMS (ACE2 monomer)",main="binding affinity");model <- lm(y~x);abline(a=0,b=1,lty=2,col="red");legend("topleft",legend=paste("R2: ",round(summary(model)$r.squared,3),sep=""),bty="n")
```

<img src="collapse_scores_files/figure-gfm/plot_correlations_v_og_DMS-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/correlations_Wuhan-Hu-1_OG-v-new-dms.pdf",sep=""),useDingbats=F))
```

## Heatmaps!

Order factor variables for plotting

``` r
#order target by order given in config
dt_final$target <- factor(dt_final$target,levels=c("Wuhan_Hu_1","E484K","N501Y","B1351","Delta"))
#order mutant as a factor for grouping by rough biochemical grouping
dt_final$mutant <- factor(dt_final$mutant, levels=c("C","P","G","V","M","L","I","A","F","W","Y","T","S","N","Q","E","D","H","K","R"))
#add character vector indicating wildtype to use as plotting symbols for wt
dt_final[,wildtype_indicator := ""]
dt_final[as.character(mutant)==as.character(wildtype),wildtype_indicator := "x"]

#make temp long-form data frame
temp <- data.table::melt(dt_final[, .(target,position,mutant,bind,delta_bind,expr,delta_expr,wildtype_indicator)],id.vars=c("target","position","mutant","wildtype_indicator"),measure.vars=c("bind","delta_bind","expr","delta_expr"),variable.name="measurement",value.name="value")

#for method to duplicate aa labels on right side of plot https://github.com/tidyverse/ggplot2/issues/3171
guide_axis_label_trans <- function(label_trans = identity, ...) {
  axis_guide <- guide_axis(...)
  axis_guide$label_trans <- rlang::as_function(label_trans)
  class(axis_guide) <- c("guide_axis_trans", class(axis_guide))
  axis_guide
}

guide_train.guide_axis_trans <- function(x, ...) {
  trained <- NextMethod()
  trained$key$.label <- x$label_trans(trained$key$.label)
  trained
}
```

Make heatmaps faceted by target, showing raw affinity and delta-affinity
of muts relative to respective

``` r
p1 <- ggplot(temp[measurement=="bind",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,11),na.value="yellow")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,12),values=c(0,1/7,7/7),na.value="yellow")+ #three notches in case I want to 'censor' closer to the 5 boundary condition
  scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,530,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_log10Ka-by-target-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_log10Ka-by-target.pdf",sep="")))
```

Second, illustrating delta_log10Ka grouped by SSM position.

``` r
p1 <- ggplot(temp[measurement=="delta_bind",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-5,2),values=c(0/7,1/7,3/7,5/7,6/7,7/7),na.value="yellow")+
  scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,530,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_delta-log10Ka-by-target-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_delta-log10Ka-by-target.pdf",sep="")))
```

Make heatmaps faceted by target, showing raw expression and
delta-expression of muts relative to respective wildtype

``` r
p1 <- ggplot(temp[measurement=="expr",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,11),na.value="yellow")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#FFFFFF","#003366"),limits=c(5,11.2),values=c(0,1/7,7/7),na.value="yellow")+ #three notches in case I want to 'censor' closer to the 5 boundary condition
  scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,530,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_expression-by-target-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_expression-by-target.pdf",sep="")))
```

Second, illustrating delta_expression grouped by SSM position.

``` r
p1 <- ggplot(temp[measurement=="delta_expr",],aes(position,mutant))+geom_tile(aes(fill=value),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#A94E35","#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-5.5,1),values=c(0/6.5,1.5/6.5,3.5/6.5,5.5/6.5,6/6.5,6.5/6.5),na.value="yellow")+
  scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,530,by=5)))+
  labs(x="",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10),axis.text.y=element_text(face="bold",size=10))+
  facet_wrap(~target,nrow=5)+
  guides(y.sec=guide_axis_label_trans())+
  geom_text(aes(label=wildtype_indicator),size=2,color="gray10")+
  theme(strip.text.x = element_text(size = 18))

p1
```

<img src="collapse_scores_files/figure-gfm/heatmap_DMS_delta-expression-by-target-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/heatmap_SSM_delta-expression-by-target.pdf",sep="")))
```

That’s the data! Other analyses in additional notebooks

Save output files.

``` r
dt_final[,.(target,wildtype,position,mutant,mutation,bind,delta_bind,n_bc_bind,n_libs_bind,bind_rep1,bind_rep2,bind_rep3,expr,delta_expr,n_bc_expr,n_libs_expr,expr_rep1,expr_rep2)] %>%
  mutate_if(is.numeric, round, digits=5) %>%
  write.csv(file=config$final_variant_scores_mut_file, row.names=F,quote=F)
```
