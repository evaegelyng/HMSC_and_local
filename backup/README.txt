Title of the study: Will be added upon acceptance of the manuscript, as the repository is public

Short summary of the study: Will be added upon acceptance of the manuscript, as the repository is public

Package versions, HPC: 

•  ape 5.7_1
•  colorspace 2.1_1
•  dplyr 1.1.3
•  fields 15.2
•  hmsc 3.0_13
•  jsonify 1.2.2
•  patchwork 1.2.0
•  phyloseq 1.50.0
•  RColorBrewer 1.1_3
•  tidyverse 2.0.0
•  vioplot 0.5.1

Package versions, local: 
•  adegenet 2.1.11
•  ape 5.8-1
•  car 3.1-3
•  circlize 0.4.16
•  ComplexHeatmap 2.22.0
•  crosstalk 1.2.1
•  dplyr 1.1.4
•  e1071 1.7-16
•  eulerr 7.0.2
•  fsa 0.10.0
•  future.apply 1.11.3
•  ggforce 0.4.2
•  ggplot2 3.5.2
•  ggpmisc 0.6.1
•  ggpubr 0.6.0
•  ggrepel 0.9.6
•  ggvenn 0.1.10
•  ggVennDiagram 1.5.2
•  gridExtra 2.3
•  hrbrthemes 0.8.7
•  htmltools 0.5.8.1
•  leaflet 2.2.2
•  leaflet.minicharts 0.6.2
•  leaflegend 1.2.1
•  lme4 1.1-37
•  magrittr 2.0.3
•  mapview 2.11.2
•  MASS 7.3-64
•  mmod 1.3.3
•  patchwork 1.3.0
•  pegas 1.3
•  phyloseq 1.50.0
•  plyr 1.8.9
•  png 0.1-8
•  RColorBrewer 1.1-3
•  readxl 1.4.5
•  scales 1.4.0
•  stringr 1.5.1
•  tibble 3.3.0
•  tidyr 1.3.1
•  tidyverse 2.0.0
•  vegan 2.6-10
•  viridis 0.6.5

Overview of folders/files and their contents: For each barcode (18S, COI, 16S), there is a folder containing the scripts that were used to run the HMSC models on a high-performance computing cluster. 
In addition, there is a folder named "Local_scripts" containing scripts for all three barcodes, which were run on a laptop.

Instructions: The scripts should be run in the following order:

1. Local_scripts/BARCODE/BARCODE_analyses.Rmd
2. BARCODE/hmsc_BARCODE_sed.r and BARCODE/hmsc_BARCODE_water.r 
3. BARCODE/Diagnostics_sed.r and BARCODE/Diagnostics_wat.r
4. BARCODE/Coef_plot_BARCODE_sed.r and BARCODE/Coef_plot_BARCODE_wat.r
5. Remaining scripts in Local_scripts*

*Within the folder Local_scripts/COI_97/ASV_level, the scripts should be run in the following order:

1. ASV_maps.Rmd
2. hap_networks.Rmd (can also be run as the last script)
3. create_data_haplot_presabs_loc_sed.R & create_data_haplot_presabs_loc_wat.R
4. ComplexHeatmap.R & Water_Sediment_NMDS_Djost.R
