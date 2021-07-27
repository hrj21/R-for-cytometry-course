# Install packages --------------------------------------------------------
install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("openCyto")
BiocManager::install("ggcyto")
devtools::install_github("jmeskas/flowCut")
devtools::install_github("DillonHammill/CytoExploreR")

library(flowCore)
library(ggcyto)

# Read individual fcs file ------------------------------------------------
my_fcs <- read.FCS("data/raw/OpA_I2_C2_IM1_WB1_R1.fcs")
my_fcs # flowframe
summary(my_fcs)
exprs(my_fcs)[1:10, ] # row per event, column per parameter
keyword(my_fcs)
names(keyword(my_fcs))
keyword(my_fcs)$`$SPILLOVER`

# Clean data using flowCut ------------------------------------------------
library(flowCut)
cut_result <- flowCut(my_fcs) #this runs flowCut and outputs diagnostics to /flowCut folder
cut_result #to see an overview
cut_result$frame # for the flowFrame to analyse
cut_result$data #for the metadata
cut_result$worstChan #for just the worst channel, useful to visalise any problems by using a plotting tool

# Clean data using flowAI -------------------------------------------------
library(flowAI)
ai_fcs <- flow_auto_qc(my_fcs)

# Store clean data in fresh flowframe -------------------------------------
cut_fcs <- cut_result$frame

# Apply instrument compensation -------------------------------------------
spillover(cut_fcs)
comp_fcs <- compensate(cut_fcs, spillover(cut_fcs)[[1]])

# Transform columns -------------------------------------------------------
to_transform <- colnames(comp_fcs)[4:9]
logicle_trans <- estimateLogicle(comp_fcs, to_transform)
trans_fcs <- transform(comp_fcs, logicle_trans)

# Plotting data -----------------------------------------------------------
ggcyto(trans_fcs, aes(`FSC-A`, `SSC-A`)) +
    geom_hex(bins = 300)

ggcyto(trans_fcs, aes(`FITC-A`, `PerCP-Cy5-5-A`)) +
    geom_hex(bins = 300)

# Working with flowsets ---------------------------------------------------
files <- list.files("data/raw", pattern = ".fcs", full.names = TRUE)
fs <- read.flowSet(files[1:4])
fs
fs[[1]]

fsApply(fs, dim)

# Cleaning data -----------------------------------------------------------
cut_fs <- fsApply(fs, function(x) {
    flow_cut_res <- flowCut(x)
    flow_cut_res$frame
})

comp_matrix <- spillover(cut_fs[[1]])[[1]]

comp_fs <- compensate(cut_fs, comp_matrix)

# Transform the flowset ---------------------------------------------------
trans_fs <- transform(comp_fs, logicle_trans)

# Plot the flowset --------------------------------------------------------
ggcyto(trans_fs, aes(`FSC-A`, `SSC-A`)) +
    geom_hex(bins = 200)

ggcyto(trans_fs, aes(`FITC-A`, `PerCP-Cy5-5-A`)) +
    geom_hex(bins = 200)


# Manual gating -----------------------------------------------------------
library(flowWorkspace)
# rectangleGate
# polygonGate
# ellipsoidGate

gs <- GatingSet(trans_fs) # gating set is the holding area for all the gates

# Scatter gate ------------------------------------------------------------
rect <- rectangleGate("FSC-A" = c(5e4, 1.2e5), "SSC-A" = c(0, 8e4), filterId = "Scatter")
rect
gs_pop_add(gs, rect, parent = "root")
gs_get_pop_paths(gs)
recompute(gs) # assign each event to the gating set

ggcyto(fs, aes(`FSC-A`, `SSC-A`)) +
    geom_hex(bins = 200) +
    geom_gate(rect)

autoplot(gs, x = "FSC-A", y = "SSC-A", "Scatter", bins = 200)

gs_pop_remove(gs, "Scatter")
recompute(gs)

# Singlet gate ------------------------------------------------------------
scatter_fs <- gs_pop_get_data(gs, "Scatter")
autoplot(scatter_fs, x = "FSC-A", y = "FSC-H", bins = 200)

poly_coords <- matrix(c(4e4, 1.2e5, 1e5, 3e4, 2.5e4, 9e4, 1.2e5, 5e4), ncol = 2)
colnames(poly_coords) <- c("FSC-A", "FSC-H")
poly <- polygonGate(poly_coords, filterId = "Singlets")
gs_pop_add(gs, poly, parent = "Scatter")
gs_get_pop_paths(gs)
recompute(gs)

autoplot(gs, x = "FSC-A", y = "FSC-H", gate = "Singlets", bins = 200)
gs_pop_remove(gs, "Singlets")

# Lymphocyte gate ---------------------------------------------------------
singlet_fs <- gs_pop_get_data(gs, "Singlets")
autoplot(singlet_fs, x = "APC-Cy7-A", y = "PE-Cy7-A", bins = 200)

cov_mat <- matrix(c(0.5, 0, 0, 0.5), ncol = 2,
                  dimnames = list(c("APC-Cy7-A", "PE-Cy7-A"), 
                                  c("APC-Cy7-A", "PE-Cy7-A")))

cov_mat

means <- c("APC-Cy7-A" = 3.5, "PE-Cy7-A" = 3.8)

ellipse <- ellipsoidGate(cov_mat, mean = means, filterId = "CD3")
gs_pop_add(gs, ellipse, parent = "Singlets")
gs_get_pop_paths(gs)
recompute(gs)

autoplot(gs, x = "APC-Cy7-A", y = "PE-Cy7-A", gate = "CD3", bins = 200)
gs_pop_remove(gs, "CD3")
recompute(gs)

# CD4 and CD8 gates -------------------------------------------------------
cd3_fs <- gs_pop_get_data(gs, "CD3")
autoplot(cd3_fs, x = "FITC-A", y = "PerCP-Cy5-5-A", bins = 200)

library(openCyto)

## CD8
cd8_auto <- fsApply(cd3_fs, function(x) {
    gate_flowclust_2d(x, 
                      xChannel = "FITC-A", 
                      yChannel = "PerCP-Cy5-5-A",
                      target = c(1, 3.5), 
                      min = c(-1, 2.5), 
                      quantile = 0.9)
})

gs_pop_add(gs, cd8_auto, parent = "CD3", name = "CD8")
recompute(gs)
autoplot(gs, x = "FITC-A", y = "PerCP-Cy5-5-A", gate = "CD8", bins = 200)
gs_pop_remove(gs, "CD8")

## CD4
cd4_auto <- fsApply(cd3_fs, function(x) {
    gate_flowclust_2d(x, xChannel = "FITC-A", yChannel = "PerCP-Cy5-5-A",
                      target = c(2.8, 0.5), 
                      max = c(4, 2),
                      K = 3,
                      quantile = 0.975)
})

gs_pop_add(gs, cd4_auto, parent = "CD3", name = "CD4")
recompute(gs)
autoplot(gs, x = "FITC-A", y = "PerCP-Cy5-5-A", gate = "CD4", bins = 200)
gs_pop_remove(gs, "CD4")

plot(gs)

autoplot(gs[[1]], bins = 200)

lapply(1:4, function(x) {
    png(paste0("plots/", x, "_gating.png"))
    print(autoplot(gs[[x]], bins = 200))
    dev.off()
})

# Get population statistics -----------------------------------------------
pop_stats <- gs_pop_get_stats(gs, type = "count")
pop_stats$percent <- gs_pop_get_stats(gs, type = "percent")$percent
gs_pop_get_stats(gs, nodes = NULL, type = pop.MFI)

medians <- gs_pop_get_stats(gs, type = pop.MFI)[, 3:8]
pop_stats <- cbind(pop_stats, medians)

write.csv(pop_stats, "Summary statistics.csv")

#gs_pop_get_stats(gs, nodes = "CD4", type = pop.MFI)

# Write out fcs files -----------------------------------------------------
cd3_fs <- gs_pop_get_data(gs, "CD3")
write.flowSet(cd3_fs, "output")

