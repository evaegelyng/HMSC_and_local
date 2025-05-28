# This repository is for joint species distribution modeling of eDNA data
# from the Coast_Sequence project.

# Before modeling, you need to run the scripts named "XXX_analyses" (with XXX being the barcode name)
# to get the input files for the HMSC models. Also run the script "Transformation_predictors.r" in the 
# Local_scripts/Across_barcodes folder to decide whether and how to transform the predictors.
# Then, for each barcode (18S, COI, 16S), run the R scripts with a name starting with "hmsc", 
# using the corresponding batch scripts. Then run the "Diagnostics" scripts to check for convergence, 
# followed by the "Coef_plot" scripts. Finally, run the "VP" script for each barcode in the Local_scripts folder.