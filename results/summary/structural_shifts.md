Analyze structural shift
================
Tyler Starr
1/26/2022

-   [Setup](#setup)
-   [Calculate pairwise C-alpha, beta, and gamma distances for each
    structure compared to
    WH1](#calculate-pairwise-c-alpha-beta-and-gamma-distances-for-each-structure-compared-to-wh1)
-   [Compare structural shift to functional
    shift](#compare-structural-shift-to-functional-shift)
-   [Repeat with aligned cryo-EM based local refined RBD:ACE2
    structures](#repeat-with-aligned-cryo-em-based-local-refined-rbdace2-structures)
    -   [Setup](#setup-1)
    -   [Calculate pairwise C-alpha, beta, and gamma distances for each
        structure compared to
        WH1](#calculate-pairwise-c-alpha-beta-and-gamma-distances-for-each-structure-compared-to-wh1-1)
    -   [Compare structural shift to functional
        shift](#compare-structural-shift-to-functional-shift-1)

This notebook analyzes structural shifts in the ACE2-bound RBD
structure, and compares to functional shifts in mutational effects.

``` r
require("knitr")
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","egg","bio3d","ggrepel")
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
if(!file.exists(config$structural_shifts_dir)){
  dir.create(file.path(config$structural_shifts_dir))
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
    ##  [1] ggrepel_0.8.1     bio3d_2.4-0       egg_0.4.5         gridExtra_2.3    
    ##  [5] forcats_0.4.0     stringr_1.4.0     dplyr_0.8.3       purrr_0.3.3      
    ##  [9] readr_1.3.1       tidyr_1.0.0       tibble_3.0.2      ggplot2_3.3.0    
    ## [13] tidyverse_1.3.0   data.table_1.12.8 yaml_2.2.0        knitr_1.26       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.0 xfun_0.11        haven_2.2.0      colorspace_1.4-1
    ##  [5] vctrs_0.3.1      generics_0.0.2   htmltools_0.4.0  rlang_0.4.7     
    ##  [9] pillar_1.4.5     glue_1.3.1       withr_2.1.2      DBI_1.1.0       
    ## [13] dbplyr_1.4.2     modelr_0.1.5     readxl_1.3.1     lifecycle_0.2.0 
    ## [17] munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0 rvest_0.3.5     
    ## [21] evaluate_0.14    parallel_3.6.2   fansi_0.4.0      broom_0.7.0     
    ## [25] Rcpp_1.0.3       scales_1.1.0     backports_1.1.5  jsonlite_1.6    
    ## [29] fs_1.3.1         hms_0.5.2        digest_0.6.23    stringi_1.4.3   
    ## [33] grid_3.6.2       cli_2.0.0        tools_3.6.2      magrittr_1.5    
    ## [37] crayon_1.3.4     pkgconfig_2.0.3  ellipsis_0.3.0   xml2_1.2.2      
    ## [41] reprex_0.3.0     lubridate_1.7.4  assertthat_0.2.1 rmarkdown_2.0   
    ## [45] httr_1.4.1       rstudioapi_0.10  R6_2.4.1         compiler_3.6.2

Define colorblind-friendly palette

``` r
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")
# The palette with black
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7")
```

## Setup

Read in PDBs. These were aligned to minimize RBD Calpha RMSD using
PyMol. For beta, omicron structures, this was performed for both RBD
units in the asymmetric unit. (These structures also included mAbs used
in the crystallization.) RBD and ACE2s were then split into separate
objects/pdbs. The pdbs in these aligned coordinates were then output.
Structures include WH1 (6m0j), N501Y (7ekf), beta (beta1=7ekg, beta2,3
not yet published), and omicron (omicron 1,2 not yet published,
omicron3=7wbp).

``` r
WH1 <- read.pdb(file=config$pdb_WH1)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
N501Y <- read.pdb(file=config$pdb_N501Y)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
beta1 <- read.pdb(file=config$pdb_beta1)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
beta2 <- read.pdb(file=config$pdb_beta2)
beta3 <- read.pdb(file=config$pdb_beta3)
delta1 <- read.pdb(file=config$pdb_delta1)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
delta2 <- read.pdb(file=config$pdb_delta2)
omicron1 <- read.pdb(file=config$pdb_omicron1)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
omicron2 <- read.pdb(file=config$pdb_omicron2)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
omicron3 <- read.pdb(file=config$pdb_omicron3)
```

## Calculate pairwise C-alpha, beta, and gamma distances for each structure compared to WH1

Note, the Cgamma calc excludes residues with no gamma-atom (alanine and
glycine) or beta-branching (and thus a gamma1 and gamma2 atom:
isoleucine, threonine, valine)

``` r
RBD_sites <- data.table(read.csv(file=config$RBD_sites))

for(voc in c("N501Y","beta1","beta2","beta3","delta1","delta2","omicron1","omicron2","omicron3")){
  RBD_sites[,paste(voc,"_v_WH1_Ca_dist",sep=""):=as.numeric(NA)]
  for(i in 1:nrow(RBD_sites)){
    site <- RBD_sites[i,site]
    WH1_calpha <- WH1$atom[WH1$atom$resno==site & WH1$atom$elety=="CA",]
    voc_calpha <- get(voc)$atom[get(voc)$atom$resno==site & get(voc)$atom$elety=="CA",]
    if(nrow(WH1_calpha)==1 & nrow(voc_calpha)==1){
      RBD_sites[i,paste(voc,"_v_WH1_Ca_dist",sep="") := sqrt((WH1_calpha$x - voc_calpha$x)^2 + (WH1_calpha$y - voc_calpha$y)^2 + (WH1_calpha$z - voc_calpha$z)^2)] 
    }else{
      #print(paste("no distance for", voc, "site",site))
    }
  }
}
```

``` r
for(voc in c("N501Y","beta1","beta2","beta3","delta1","delta2","omicron1","omicron2","omicron3")){
  RBD_sites[,paste(voc,"_v_WH1_Cb_dist",sep=""):=as.numeric(NA)]
  for(i in 1:nrow(RBD_sites)){
    site <- RBD_sites[i,site]
    WH1_cbeta <- WH1$atom[WH1$atom$resno==site & WH1$atom$elety=="CB",]
    voc_cbeta <- get(voc)$atom[get(voc)$atom$resno==site & get(voc)$atom$elety=="CB",]
    if(nrow(WH1_cbeta)==1 & nrow(voc_cbeta)==1){
      RBD_sites[i,paste(voc,"_v_WH1_Cb_dist",sep="") := sqrt((WH1_cbeta$x - voc_cbeta$x)^2 + (WH1_cbeta$y - voc_cbeta$y)^2 + (WH1_cbeta$z - voc_cbeta$z)^2)] 
    }else{
      #print(paste("no distance for", voc, "site",site))
    }
  }
}
```

``` r
for(voc in c("N501Y","beta1","beta2","beta3","delta1","delta2","omicron1","omicron2","omicron3")){
  RBD_sites[,paste(voc,"_v_WH1_Cg_dist",sep=""):=as.numeric(NA)]
  for(i in 1:nrow(RBD_sites)){
    site <- RBD_sites[i,site]
    WH1_cgamma <- WH1$atom[WH1$atom$resno==site & WH1$atom$elety %in% c("CG","SG","OG"),]
    voc_cgamma <- get(voc)$atom[get(voc)$atom$resno==site & get(voc)$atom$elety %in% c("CG","SG","OG"),]
    if(nrow(WH1_cgamma)==1 & nrow(voc_cgamma)==1){
      RBD_sites[i,paste(voc,"_v_WH1_Cg_dist",sep="") := sqrt((WH1_cgamma$x - voc_cgamma$x)^2 + (WH1_cgamma$y - voc_cgamma$y)^2 + (WH1_cgamma$z - voc_cgamma$z)^2)] 
    }else{
      #print(paste("no distance for", voc, "site",site))
    }
  }
}
```

Average across all atoms in each residue

``` r
for(voc in c("N501Y","beta1","beta2","beta3","delta1","delta2","omicron1","omicron2","omicron3")){
  RBD_sites[,paste(voc,"_v_WH1_all_atom_dist",sep=""):=as.numeric(NA)]
  for(i in 1:nrow(RBD_sites)){
    site <- RBD_sites[i,site]
    WH1_atoms <- WH1$atom[WH1$atom$resno==site,]
    voc_atoms <- get(voc)$atom[get(voc)$atom$resno==site,]
    if(nrow(WH1_atoms)>1 & nrow(voc_atoms)>1 & nrow(WH1_atoms)==nrow(voc_atoms)){
      if(WH1_atoms$resid[1]==voc_atoms$resid[1]){
        RBD_sites[i,paste(voc,"_v_WH1_all_atom_dist",sep="") := mean(sqrt((WH1_atoms$x - voc_atoms$x)^2 + (WH1_atoms$y - voc_atoms$y)^2 + (WH1_atoms$z - voc_atoms$z)^2))] 
      }
    }
  }
}
```

Line plots of Calpha, beta, gamma, atom-averaged displacement compared
to WH1 structure

``` r
#define focal bg for others to compare to
calpha <- melt(RBD_sites[,.(site,N501Y_v_WH1_Ca_dist,beta1_v_WH1_Ca_dist,beta2_v_WH1_Ca_dist,beta3_v_WH1_Ca_dist,delta1_v_WH1_Ca_dist,delta2_v_WH1_Ca_dist)], id.vars=c("site"))
calpha[,rep:=1]
calpha[variable=="N501Y_v_WH1_Ca_dist",variable:="N501Y"]
calpha[variable=="beta1_v_WH1_Ca_dist",variable:="Beta"]
calpha[variable=="beta2_v_WH1_Ca_dist",rep:=2]
calpha[variable=="beta2_v_WH1_Ca_dist",variable:="Beta"]
calpha[variable=="beta3_v_WH1_Ca_dist",rep:=3]
calpha[variable=="beta3_v_WH1_Ca_dist",variable:="Beta"]
calpha[variable=="delta1_v_WH1_Ca_dist",variable:="Delta"]
calpha[variable=="delta2_v_WH1_Ca_dist",rep:=2]
calpha[variable=="delta2_v_WH1_Ca_dist",variable:="Delta"]

#define colors for each bg
group.colors <- c("Wuhan-Hu-1" = cbPalette[1], "N501Y" = cbPalette[3], "E484K" = cbPalette[5], "Beta"=cbPalette[6], "Delta"=cbPalette[7], "Omicron"=cbPalette[2])

#define order for plotting of bgs
calpha$variable <- factor(calpha$variable,levels=c("E484K","N501Y","Beta","Delta"))

ggplot(data=calpha[rep==1 & site %in% seq(334,515),], aes(x=site, y=value, color=variable))+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("C-alpha distance versus Wuhan-Hu-1 [6m0j] (A)")+
  geom_text_repel(aes(label=ifelse(((value > 0.75) & variable=="Beta"),as.character(site),'')),size=3,color="gray40")
```

<img src="structural_shifts_files/figure-gfm/line_plots_Calpha-dist_v_WH1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/Calpha_v_WH1.pdf",sep=""),useDingbats=F))
```

``` r
#define focal bg for others to compare to
cbeta <- melt(RBD_sites[,.(site,N501Y_v_WH1_Cb_dist,beta1_v_WH1_Cb_dist,beta2_v_WH1_Cb_dist,beta3_v_WH1_Cb_dist,delta1_v_WH1_Cb_dist,delta2_v_WH1_Cb_dist)], id.vars=c("site"))
cbeta[,rep:=1]
cbeta[variable=="N501Y_v_WH1_Cb_dist",variable:="N501Y"]
cbeta[variable=="beta1_v_WH1_Cb_dist",variable:="Beta"]
cbeta[variable=="beta2_v_WH1_Cb_dist",rep:=2]
cbeta[variable=="beta2_v_WH1_Cb_dist",variable:="Beta"]
cbeta[variable=="beta3_v_WH1_Cb_dist",rep:=3]
cbeta[variable=="beta3_v_WH1_Cb_dist",variable:="Beta"]
cbeta[variable=="delta1_v_WH1_Cb_dist",variable:="Delta"]
cbeta[variable=="delta2_v_WH1_Cb_dist",rep:=2]
cbeta[variable=="delta2_v_WH1_Cb_dist",variable:="Delta"]

#define order for plotting of bgs
cbeta$variable <- factor(cbeta$variable,levels=c("E484K","N501Y","Beta","Delta"))

ggplot(data=cbeta[rep==1 & site %in% seq(334,515),], aes(x=site, y=value, color=variable))+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("C-beta distance versus Wuhan-Hu-1 [6m0j] (A)")
```

<img src="structural_shifts_files/figure-gfm/line_plots_Cbeta-dist_v_WH1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/Cbeta_v_WH1.pdf",sep=""),useDingbats=F))
```

``` r
#define focal bg for others to compare to
cgamma <- melt(RBD_sites[,.(site,N501Y_v_WH1_Cg_dist,beta1_v_WH1_Cg_dist,beta2_v_WH1_Cg_dist,beta3_v_WH1_Cg_dist,delta1_v_WH1_Cg_dist,delta2_v_WH1_Cg_dist)], id.vars=c("site"))
cgamma[,rep:=1]
cgamma[variable=="N501Y_v_WH1_Cg_dist",variable:="N501Y"]
cgamma[variable=="beta1_v_WH1_Cg_dist",variable:="Beta"]
cgamma[variable=="beta2_v_WH1_Cg_dist",rep:=2]
cgamma[variable=="beta2_v_WH1_Cg_dist",variable:="Beta"]
cgamma[variable=="beta3_v_WH1_Cg_dist",rep:=3]
cgamma[variable=="beta3_v_WH1_Cg_dist",variable:="Beta"]
cgamma[variable=="delta1_v_WH1_Cg_dist",variable:="Delta"]
cgamma[variable=="delta2_v_WH1_Cg_dist",rep:=2]
cgamma[variable=="delta2_v_WH1_Cg_dist",variable:="Delta"]

#define order for plotting of bgs
cgamma$variable <- factor(cgamma$variable,levels=c("E484K","N501Y","Beta","Delta"))

ggplot(data=cgamma[rep==1 & site %in% seq(334,515),], aes(x=site, y=value, color=variable))+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("C-gamma distance versus Wuhan-Hu-1 [6m0j] (A)")
```

<img src="structural_shifts_files/figure-gfm/line_plots_Cgamma-dist_v_WH1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/Cgamma_v_WH1.pdf",sep=""),useDingbats=F))
```

``` r
#define focal bg for others to compare to
all_atom <- melt(RBD_sites[,.(site,N501Y_v_WH1_all_atom_dist,beta1_v_WH1_all_atom_dist,beta2_v_WH1_all_atom_dist,beta3_v_WH1_all_atom_dist,delta1_v_WH1_all_atom_dist,delta2_v_WH1_all_atom_dist)], id.vars=c("site"))
all_atom[,rep:=1]
all_atom[variable=="N501Y_v_WH1_all_atom_dist",variable:="N501Y"]
all_atom[variable=="beta1_v_WH1_all_atom_dist",variable:="Beta"]
all_atom[variable=="beta2_v_WH1_all_atom_dist",rep:=2]
all_atom[variable=="beta2_v_WH1_all_atom_dist",variable:="Beta"]
all_atom[variable=="beta3_v_WH1_all_atom_dist",rep:=3]
all_atom[variable=="beta3_v_WH1_all_atom_dist",variable:="Beta"]
all_atom[variable=="delta1_v_WH1_all_atom_dist",variable:="Delta"]
all_atom[variable=="delta2_v_WH1_all_atom_dist",rep:=2]
all_atom[variable=="delta2_v_WH1_all_atom_dist",variable:="Delta"]

#define order for plotting of bgs
all_atom$variable <- factor(all_atom$variable,levels=c("E484K","N501Y","Beta","Delta"))

ggplot(data=all_atom[rep==1 & site %in% seq(334,515),], aes(x=site, y=value, color=variable))+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("Average atomic displacement versus Wuhan-Hu-1 [6m0j] (A)")+
  geom_text_repel(aes(label=ifelse(((value > 1) & variable=="Beta"),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 6 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/line_plots_average-atom-dist_v_WH1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/all_atom_v_WH1.pdf",sep=""),useDingbats=F))
```

## Compare structural shift to functional shift

We have a metric of functional perturbation at a site, derived from our
deep mutational scanning measurements. Let???s see how structural
perturbations correlate (or don???t) with structural perturbations.

``` r
JSD_bind <- data.table(read.csv(file=config$JSD_v_WH1_file, stringsAsFactors=F))
JSD_expr <- data.table(read.csv(file=config$JSD_v_WH1_expr_file, stringsAsFactors=F))

for(i in 1:nrow(calpha)){
  calpha$JSD_bind[i] <- JSD_bind[site==calpha$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==calpha$variable[i],JSD_min3bc]
  calpha$JSD_expr[i] <- JSD_expr[site==calpha$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==calpha$variable[i],JSD_min3bc]
}

for(i in 1:nrow(cbeta)){
  cbeta$JSD_bind[i] <- JSD_bind[site==cbeta$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==cbeta$variable[i],JSD_min3bc]
  cbeta$JSD_expr[i] <- JSD_expr[site==cbeta$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==cbeta$variable[i],JSD_min3bc]
}

for(i in 1:nrow(cgamma)){
  cgamma$JSD_bind[i] <- JSD_bind[site==cgamma$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==cgamma$variable[i],JSD_min3bc]
  cgamma$JSD_expr[i] <- JSD_expr[site==cgamma$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==cgamma$variable[i],JSD_min3bc]
}

for(i in 1:nrow(all_atom)){
  all_atom$JSD_bind[i] <- JSD_bind[site==all_atom$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==all_atom$variable[i],JSD_min3bc]
  all_atom$JSD_expr[i] <- JSD_expr[site==all_atom$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==all_atom$variable[i],JSD_min3bc]
}
```

Calpha displacement:

``` r
ggplot(data=calpha[rep==1 & site %in% seq(334,515),],aes(x=abs(value),y=JSD_bind,color=variable))+
  geom_point(pch=19)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  xlab("C-alpha displacement versus Wuhan-Hu-1 [6m0j]")+
  ylab("JS divergence versus Wuhan-Hu-1")+
  facet_wrap(~variable,nrow=1)+
  theme(strip.text.x = element_text(size = 18))+
  geom_text_repel(aes(label=ifelse(((JSD_bind > 0.1 & abs(value) > 0.75) | (JSD_bind > 0.15 & abs(value) > 0.25)),as.character(site),'')),size=3,color="gray40")
```

<img src="structural_shifts_files/figure-gfm/scatter_plots_calpha_v_JSD_v_WH1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/scatter_JSD_v_calpha.pdf",sep=""),useDingbats=F))
```

Cbeta displacement:

``` r
ggplot(data=cbeta[rep==1 & site %in% seq(334,515),],aes(x=abs(value),y=JSD_bind,color=variable))+
  geom_point(pch=19)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  xlab("C-beta displacement versus Wuhan-Hu-1 [6m0j]")+
  ylab("JS divergence versus Wuhan-Hu-1")+
  facet_wrap(~variable,nrow=1)+
  theme(strip.text.x = element_text(size = 18))+
  geom_text_repel(aes(label=ifelse(((JSD_bind > 0.1 & abs(value) > 0.75) | (JSD_bind > 0.15 & abs(value) > 0.25)),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 42 rows containing missing values (geom_point).

    ## Warning: Removed 42 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/scatter_plots_cbeta_v_JSD_v_WH1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/scatter_JSD_v_cbeta.pdf",sep=""),useDingbats=F))
```

C-gamma displacement:

``` r
ggplot(data=cgamma[rep==1 & site %in% seq(334,515),],aes(x=abs(value),y=JSD_bind,color=variable))+
  geom_point(pch=19)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  xlab("C-gamma displacement versus Wuhan-Hu-1 [6m0j]")+
  ylab("JS divergence versus Wuhan-Hu-1")+
  facet_wrap(~variable,nrow=1)+
  theme(strip.text.x = element_text(size = 18))+
  geom_text_repel(aes(label=ifelse(((JSD_bind > 0.1 & abs(value) > 0.75) | (JSD_bind > 0.15 & abs(value) > 0.25)),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 165 rows containing missing values (geom_point).

    ## Warning: Removed 165 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/scatter_plots_cgamma_v_JSD_v_WH1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/scatter_JSD_v_cgamma.pdf",sep=""),useDingbats=F))
```

All-atom displacement:

``` r
ggplot(data=all_atom[rep==1 & site %in% seq(334,515),],aes(x=abs(value),y=JSD_bind,color=variable))+
  geom_point(pch=19)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  xlab("Average atomic displacement versus Wuhan-Hu-1 [6m0j]")+
  ylab("JS divergence versus Wuhan-Hu-1")+
  facet_wrap(~variable,nrow=1)+
  theme(strip.text.x = element_text(size = 18))+
  geom_text_repel(aes(label=ifelse(((JSD_bind > 0.1 & abs(value) > 0.75) | (JSD_bind > 0.15 & abs(value) > 0.25)),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 6 rows containing missing values (geom_point).

    ## Warning: Removed 6 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/scatter_plots_all-atom_v_JSD_v_WH1-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/scatter_JSD_v_all_atom.pdf",sep=""),useDingbats=F))
```

# Repeat with aligned cryo-EM based local refined RBD:ACE2 structures

## Setup

Read in PDBs. These were aligned to minimize RBD Calpha RMSD using
PyMol. Structures include WH1 (7kmb), N501Y (7ekf), beta (7vx4), and
delta (7v8b)

``` r
WH1 <- read.pdb(file=config$pdb_WH1_cryo)
alpha <- read.pdb(file=config$pdb_alpha_cryo)
beta <- read.pdb(file=config$pdb_beta_cryo)
```

    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
delta <- read.pdb(file=config$pdb_delta_cryo)
```

## Calculate pairwise C-alpha, beta, and gamma distances for each structure compared to WH1

Note, the Cgamma calc excludes residues with no gamma-atom (alanine and
glycine) or beta-branching (and thus a gamma1 and gamma2 atom:
isoleucine, threonine, valine)

``` r
RBD_sites <- data.table(read.csv(file=config$RBD_sites))

for(voc in c("alpha","beta","delta")){
  RBD_sites[,paste(voc,"_v_WH1_Ca_dist",sep=""):=as.numeric(NA)]
  for(i in 1:nrow(RBD_sites)){
    site <- RBD_sites[i,site]
    WH1_calpha <- WH1$atom[WH1$atom$resno==site & WH1$atom$elety=="CA",]
    voc_calpha <- get(voc)$atom[get(voc)$atom$resno==site & get(voc)$atom$elety=="CA",]
    if(nrow(WH1_calpha)==1 & nrow(voc_calpha)==1){
      RBD_sites[i,paste(voc,"_v_WH1_Ca_dist",sep="") := sqrt((WH1_calpha$x - voc_calpha$x)^2 + (WH1_calpha$y - voc_calpha$y)^2 + (WH1_calpha$z - voc_calpha$z)^2)] 
    }else{
      #print(paste("no distance for", voc, "site",site))
    }
  }
}
```

``` r
for(voc in c("alpha","beta","delta")){
  RBD_sites[,paste(voc,"_v_WH1_Cb_dist",sep=""):=as.numeric(NA)]
  for(i in 1:nrow(RBD_sites)){
    site <- RBD_sites[i,site]
    WH1_cbeta <- WH1$atom[WH1$atom$resno==site & WH1$atom$elety=="CB",]
    voc_cbeta <- get(voc)$atom[get(voc)$atom$resno==site & get(voc)$atom$elety=="CB",]
    if(nrow(WH1_cbeta)==1 & nrow(voc_cbeta)==1){
      RBD_sites[i,paste(voc,"_v_WH1_Cb_dist",sep="") := sqrt((WH1_cbeta$x - voc_cbeta$x)^2 + (WH1_cbeta$y - voc_cbeta$y)^2 + (WH1_cbeta$z - voc_cbeta$z)^2)] 
    }else{
      #print(paste("no distance for", voc, "site",site))
    }
  }
}
```

``` r
for(voc in c("alpha","beta","delta")){
  RBD_sites[,paste(voc,"_v_WH1_Cg_dist",sep=""):=as.numeric(NA)]
  for(i in 1:nrow(RBD_sites)){
    site <- RBD_sites[i,site]
    WH1_cgamma <- WH1$atom[WH1$atom$resno==site & WH1$atom$elety %in% c("CG","SG","OG"),]
    voc_cgamma <- get(voc)$atom[get(voc)$atom$resno==site & get(voc)$atom$elety %in% c("CG","SG","OG"),]
    if(nrow(WH1_cgamma)==1 & nrow(voc_cgamma)==1){
      RBD_sites[i,paste(voc,"_v_WH1_Cg_dist",sep="") := sqrt((WH1_cgamma$x - voc_cgamma$x)^2 + (WH1_cgamma$y - voc_cgamma$y)^2 + (WH1_cgamma$z - voc_cgamma$z)^2)] 
    }else{
      #print(paste("no distance for", voc, "site",site))
    }
  }
}
```

Average across all atoms in each residue

``` r
for(voc in c("alpha","beta","delta")){
  RBD_sites[,paste(voc,"_v_WH1_all_atom_dist",sep=""):=as.numeric(NA)]
  for(i in 1:nrow(RBD_sites)){
    site <- RBD_sites[i,site]
    WH1_atoms <- WH1$atom[WH1$atom$resno==site,]
    voc_atoms <- get(voc)$atom[get(voc)$atom$resno==site,]
    if(nrow(WH1_atoms)>1 & nrow(voc_atoms)>1 & nrow(WH1_atoms)==nrow(voc_atoms)){
      if(WH1_atoms$resid[1]==voc_atoms$resid[1]){
        RBD_sites[i,paste(voc,"_v_WH1_all_atom_dist",sep="") := mean(sqrt((WH1_atoms$x - voc_atoms$x)^2 + (WH1_atoms$y - voc_atoms$y)^2 + (WH1_atoms$z - voc_atoms$z)^2))] 
      }
    }
  }
}
```

Line plots of Calpha, beta, gamma displacement compared to WH1 structure

``` r
#define focal bg for others to compare to
calpha <- melt(RBD_sites[,.(site,alpha_v_WH1_Ca_dist,beta_v_WH1_Ca_dist,delta_v_WH1_Ca_dist)], id.vars=c("site"))
calpha[variable=="alpha_v_WH1_Ca_dist",variable:="N501Y"]
calpha[variable=="beta_v_WH1_Ca_dist",variable:="Beta"]
calpha[variable=="delta_v_WH1_Ca_dist",variable:="Delta"]

#define colors for each bg
group.colors <- c("Wuhan-Hu-1" = cbPalette[1], "N501Y" = cbPalette[3], "E484K" = cbPalette[5], "Beta"=cbPalette[6], "Delta"=cbPalette[7], "Omicron"=cbPalette[2])

#define order for plotting of bgs
calpha$variable <- factor(calpha$variable,levels=c("E484K","N501Y","Beta","Delta"))

ggplot(data=calpha[site %in% seq(334,515),], aes(x=site, y=value, color=variable))+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("C-alpha distance versus Wuhan-Hu-1 [7kmb] (A)")+
  geom_text_repel(aes(label=ifelse(((value > 1) & variable=="Beta"),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 3 row(s) containing missing values (geom_path).

    ## Warning: Removed 3 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/line_plots_Calpha-dist_v_WH1_cryo-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/Calpha_v_WH1_cryo.pdf",sep=""),useDingbats=F))
```

``` r
#define focal bg for others to compare to
cbeta <- melt(RBD_sites[,.(site,alpha_v_WH1_Cb_dist,beta_v_WH1_Cb_dist,delta_v_WH1_Cb_dist)], id.vars=c("site"))
cbeta[variable=="alpha_v_WH1_Cb_dist",variable:="N501Y"]
cbeta[variable=="beta_v_WH1_Cb_dist",variable:="Beta"]
cbeta[variable=="delta_v_WH1_Cb_dist",variable:="Delta"]

#define order for plotting of bgs
cbeta$variable <- factor(cbeta$variable,levels=c("E484K","N501Y","Beta","Delta"))

ggplot(data=cbeta[site %in% seq(334,515),], aes(x=site, y=value, color=variable))+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("C-beta distance versus Wuhan-Hu-1 [7kmb] (A)")+
  geom_text_repel(aes(label=ifelse(((value > 2) & variable=="Beta"),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 3 row(s) containing missing values (geom_path).

    ## Warning: Removed 45 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/line_plots_Cbeta-dist_v_WH1_cryo-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/Cbeta_v_WH1_cryo.pdf",sep=""),useDingbats=F))
```

``` r
#define focal bg for others to compare to
cgamma <- melt(RBD_sites[,.(site,alpha_v_WH1_Cg_dist,beta_v_WH1_Cg_dist,delta_v_WH1_Cg_dist)], id.vars=c("site"))
cgamma[variable=="alpha_v_WH1_Cg_dist",variable:="N501Y"]
cgamma[variable=="beta_v_WH1_Cg_dist",variable:="Beta"]
cgamma[variable=="delta_v_WH1_Cg_dist",variable:="Delta"]

#define order for plotting of bgs
cgamma$variable <- factor(cgamma$variable,levels=c("E484K","N501Y","Beta","Delta"))

ggplot(data=cgamma[site %in% seq(334,515),], aes(x=site, y=value, color=variable))+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("C-gamma distance versus Wuhan-Hu-1 [7kmb] (A)")
```

    ## Warning: Removed 3 row(s) containing missing values (geom_path).

<img src="structural_shifts_files/figure-gfm/line_plots_Cgamma-dist_v_WH1_cryo-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/Cgamma_v_WH1_cryo.pdf",sep=""),useDingbats=F))
```

``` r
#define focal bg for others to compare to
all_atom <- melt(RBD_sites[,.(site,alpha_v_WH1_all_atom_dist,beta_v_WH1_all_atom_dist,delta_v_WH1_all_atom_dist)], id.vars=c("site"))
all_atom[variable=="alpha_v_WH1_all_atom_dist",variable:="N501Y"]
all_atom[variable=="beta_v_WH1_all_atom_dist",variable:="Beta"]
all_atom[variable=="delta_v_WH1_all_atom_dist",variable:="Delta"]

#define order for plotting of bgs
all_atom$variable <- factor(all_atom$variable,levels=c("E484K","N501Y","Beta","Delta"))

ggplot(data=all_atom[site %in% seq(334,515),], aes(x=site, y=value, color=variable))+
  geom_line(size=1)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  scale_x_continuous(expand=c(0.01,0.01),breaks=c(331,seq(335,530,by=5)))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold",size=10))+
  ylab("Average atomic displacement versus Wuhan-Hu-1 [7kmb] (A)")+
  geom_text_repel(aes(label=ifelse(((value > 2) & variable=="Beta"),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 3 row(s) containing missing values (geom_path).

    ## Warning: Removed 9 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/line_plots_all_atom-dist_v_WH1_cryo-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/all_atom_v_WH1_cryo.pdf",sep=""),useDingbats=F))
```

## Compare structural shift to functional shift

We have a metric of functional perturbation at a site, derived from our
deep mutational scanning measurements. Let???s see how structural
perturbations correlate (or don???t) with structural perturbations.

``` r
JSD_bind <- data.table(read.csv(file=config$JSD_v_WH1_file, stringsAsFactors=F))
JSD_expr <- data.table(read.csv(file=config$JSD_v_WH1_expr_file, stringsAsFactors=F))

for(i in 1:nrow(calpha)){
  calpha$JSD_bind[i] <- JSD_bind[site==calpha$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==calpha$variable[i],JSD_min3bc]
  calpha$JSD_expr[i] <- JSD_expr[site==calpha$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==calpha$variable[i],JSD_min3bc]
}

for(i in 1:nrow(cbeta)){
  cbeta$JSD_bind[i] <- JSD_bind[site==cbeta$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==cbeta$variable[i],JSD_min3bc]
  cbeta$JSD_expr[i] <- JSD_expr[site==cbeta$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==cbeta$variable[i],JSD_min3bc]
}

for(i in 1:nrow(cgamma)){
  cgamma$JSD_bind[i] <- JSD_bind[site==cgamma$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==cgamma$variable[i],JSD_min3bc]
  cgamma$JSD_expr[i] <- JSD_expr[site==cgamma$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==cgamma$variable[i],JSD_min3bc]
}

for(i in 1:nrow(all_atom)){
  all_atom$JSD_bind[i] <- JSD_bind[site==all_atom$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==all_atom$variable[i],JSD_min3bc]
  all_atom$JSD_expr[i] <- JSD_expr[site==all_atom$site[i] & bg_1=="Wuhan-Hu-1" & bg_2==all_atom$variable[i],JSD_min3bc]
}
```

Calpha displacement:

``` r
ggplot(data=calpha[site %in% seq(334,515),],aes(x=abs(value),y=JSD_bind,color=variable))+
  geom_point(pch=19)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  xlab("C-alpha displacement versus Wuhan-Hu-1 [7kmb]")+
  ylab("JS divergence versus Wuhan-Hu-1")+
  facet_wrap(~variable,nrow=1)+
  theme(strip.text.x = element_text(size = 18))+
  geom_text_repel(aes(label=ifelse(((JSD_bind > 0.1 & abs(value) > 0.75) | (JSD_bind > 0.15 & abs(value) > 0.25)),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 3 rows containing missing values (geom_point).

    ## Warning: Removed 3 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/scatter_plots_calpha_v_JSD_v_WH1_cryo-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/scatter_JSD_v_calpha_cryo.pdf",sep=""),useDingbats=F))
```

Cbeta displacement:

``` r
ggplot(data=cbeta[site %in% seq(334,515),],aes(x=abs(value),y=JSD_bind,color=variable))+
  geom_point(pch=19)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  xlab("C-beta displacement versus Wuhan-Hu-1 [7kmb]")+
  ylab("JS divergence versus Wuhan-Hu-1")+
  facet_wrap(~variable,nrow=1)+
  theme(strip.text.x = element_text(size = 18))+
  geom_text_repel(aes(label=ifelse(((JSD_bind > 0.1 & abs(value) > 0.75) | (JSD_bind > 0.15 & abs(value) > 0.25)),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 45 rows containing missing values (geom_point).

    ## Warning: Removed 45 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/scatter_plots_cbeta_v_JSD_v_WH1_cryo-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/scatter_JSD_v_cbeta_cryo.pdf",sep=""),useDingbats=F))
```

C-gamma displacement:

``` r
ggplot(data=cgamma[site %in% seq(334,515),],aes(x=abs(value),y=JSD_bind,color=variable))+
  geom_point(pch=19)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  xlab("C-gamma displacement versus Wuhan-Hu-1 [7kmb]")+
  ylab("JS divergence versus Wuhan-Hu-1")+
  facet_wrap(~variable,nrow=1)+
  theme(strip.text.x = element_text(size = 18))+
  geom_text_repel(aes(label=ifelse(((JSD_bind > 0.1 & abs(value) > 0.75) | (JSD_bind > 0.15 & abs(value) > 0.25)),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 168 rows containing missing values (geom_point).

    ## Warning: Removed 168 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/scatter_plots_cgamma_v_JSD_v_WH1_cryo-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/scatter_JSD_v_cgamma_cryo.pdf",sep=""),useDingbats=F))
```

``` r
ggplot(data=all_atom[site %in% seq(334,515),],aes(x=abs(value),y=JSD_bind,color=variable))+
  geom_point(pch=19)+
  scale_color_manual(values=group.colors)+
  theme_classic()+
  xlab("Average atomic displacement versus Wuhan-Hu-1 [7kmb]")+
  ylab("JS divergence versus Wuhan-Hu-1")+
  facet_wrap(~variable,nrow=1)+
  theme(strip.text.x = element_text(size = 18))+
  geom_text_repel(aes(label=ifelse(((JSD_bind > 0.1 & abs(value) > 0.75) | (JSD_bind > 0.15 & abs(value) > 0.25)),as.character(site),'')),size=3,color="gray40")
```

    ## Warning: Removed 9 rows containing missing values (geom_point).

    ## Warning: Removed 9 rows containing missing values (geom_text_repel).

<img src="structural_shifts_files/figure-gfm/scatter_plots_all_atom_v_JSD_v_WH1_cryo-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$structural_shifts_dir,"/scatter_JSD_v_all_atom_cryo.pdf",sep=""),useDingbats=F))
```
