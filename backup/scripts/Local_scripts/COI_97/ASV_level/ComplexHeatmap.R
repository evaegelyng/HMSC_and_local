# libraries
library(vegan)
library(ComplexHeatmap)
library(circlize)

# load water data
nmds_water <- read.table("../../../Tekstfiler/COI/COI_ASV/wat_Djost.txt", sep=",", header=TRUE)
row.names(nmds_water) <- nmds_water$X
nmds_water <- nmds_water[,-1] # remove "X" column
nmds_water[is.na(nmds_water)] = 0 # set NA values to 0
#diag(nmds_water)=NA # use NA for self-comparisons

# Now we load the water data on number of comparisons
mat_values_water = read.table("../../../Tekstfiler/COI/COI_ASV/water_no_of_comp.txt", sep="\t", header=TRUE)
mat_values_water[is.na(mat_values_water)] = 0 # set NA values to 0

# Color scheme
col_fun_water = colorRamp2(c(0, 0.46, 0.92), c("blue", "white", "red"))
col_fun_water(seq(-3, 3))

#make water plot
heatmap_wat<-grid.grabExpr(draw(Heatmap(as.matrix(nmds_water), 
        #        clustering_method_rows = "complete", # Consider average or complete (default)
        na_col="black",
        #        column_km = 7,
        #        row_km = 7,
        name="Jost's D", 
        col = col_fun_water,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", mat_values_water[i, j]), x, y, gp = gpar(fontsize = 11, fontface = 'bold'))
        }),show_heatmap_legend=TRUE))

saveRDS(heatmap_wat,"../../../Plots/Heatmap_Djost_wat.rds")

# load sediment data
nmds_sediment <- read.table("../../../Tekstfiler/COI/COI_ASV/Sed_Djost.txt", sep=",", header=TRUE)
row.names(nmds_sediment) <- nmds_sediment$X
nmds_sediment <- nmds_sediment[,-1] # remove "X" column
nmds_sediment[is.na(nmds_sediment)] = 0 # set NA values to 0
#diag(nmds_sediment)=NA # use NA for self-comparisons

# Now we load the sediment data on number of comparisons
mat_values_sediment = read.table("../../../Tekstfiler/COI/COI_ASV/Sediment_no_of_comp.txt", sep="\t", header=TRUE)
mat_values_sediment[is.na(mat_values_sediment)] = 0 # set NA values to 0

# Color scheme
col_fun_sediment = colorRamp2(c(0, 0.47, 0.94), c("blue", "white", "red"))
col_fun_sediment(seq(-3, 3))

#make Sediment plot
heatmap_sed <- grid.grabExpr(draw(Heatmap(as.matrix(nmds_sediment), 
                clustering_method_rows = "complete", # Consider average or complete (default)
        na_col="black",
        #        column_km = 7,
        #        row_km = 7,
        name="Jost's D", 
        col = col_fun_sediment,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", mat_values_sediment[i, j]), x, y, gp = gpar(fontsize = 11, fontface = 'bold'))
        }),show_heatmap_legend=FALSE))

saveRDS(heatmap_sed,"../../../Plots/Heatmap_Djost_sed.rds")

# Things to consider adding
# row_split = data.frame(rep(c("A", "B"), 31), rep(c("C", "D"), each = 15)), --> didn't make this work
