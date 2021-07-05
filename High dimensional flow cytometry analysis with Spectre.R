# Install packages --------------------------------------------------------
install.packages(c("devtools", "BiocManager"))
library(devtools)
install_github("immunedynamics/spectre")
library(BiocManager)
install("EnhancedVolcano")
library(Spectre)

## Check if all required packages have been installed
package.check()

## Load all required packages
package.load()

# Read data ---------------------------------------------------------------
list.files("data/clean", pattern = "fcs")

data_list <- read.files(file.loc = "data/clean", 
                        file.type = ".fcs", 
                        do.embed.file.names = TRUE)

# Check data --------------------------------------------------------------
do.list.summary(data_list)

# Merge data --------------------------------------------------------------
cell_dat <- do.merge.files(dat = data_list)

# Transform channel values ------------------------------------------------
names(cell_dat)
to_asinh <- names(cell_dat)[1:8]
to_asinh

cell_dat <- do.asinh(cell_dat, to_asinh, cofactor = 500)
# cell_dat <- do.logicle(cell_dat, to_asinh)
transformed_cols <- paste0(to_asinh, "_asinh")

# Plot channels -----------------------------------------------------------
lapply(transformed_cols, function(x) {
    make.colour.plot(do.subsample(cell_dat, 1e4), 
                     x, 
                     "AF700_CD45_asinh",
                     path = "plots",
                     save.to.disk = TRUE)
})

# Read metadata -----------------------------------------------------------
meta_dat <- read.csv("sample_details.csv")
meta_dat

cell_dat <- do.add.cols(cell_dat, 
                        base.col = "FileName", 
                        add.dat = meta_dat, 
                        add.by = "FileName")

cell_dat

# Batch alignment ---------------------------------------------------------
reference_files <- c("Mock_01_A", "Mock_05_B")
reference_dat <- do.filter(cell_dat, "Sample", reference_files)

cytonorm <- prep.cytonorm(dat = reference_dat,
                          cellular.cols = transformed_cols,
                          cluster.cols = transformed_cols,
                          batch.col = "Batch",
                          sample.col = "Sample") 

cytonorm_sub <- do.subsample(cytonorm$dt, 1e4)
cytonorm_sub <- run.umap(cytonorm_sub, use.cols = transformed_cols)

make.colour.plot(cytonorm_sub, 
                 'UMAP_X', 
                 'UMAP_Y', 
                 col.axis = 'Sample', 
                 col.type = 'factor', 
                 path = "plots")

make.colour.plot(cytonorm_sub, 
                 'UMAP_X', 
                 'UMAP_Y', 
                 col.axis = 'prep.fsom.metacluster', 
                 col.type = 'factor', 
                 path = "plots")

cytonorm <- train.cytonorm(model = cytonorm, align.cols = transformed_cols)

cell_dat <- run.cytonorm(dat = cell_dat, model = cytonorm, batch.col = "Batch")

aligned_cols <- paste0(transformed_cols, '_aligned')

reference_dat <- do.filter(cell_dat, "Sample", reference_files)
reference_sub <- do.subsample(reference_dat, 1e4)
reference_sub <- run.umap(reference_sub, use.cols = aligned_cols)

make.colour.plot(reference_sub, 'UMAP_X', 'UMAP_Y', "Batch", 'factor', path = "plots")

# Clustering --------------------------------------------------------------
aligned_cols <- paste0(transformed_cols, '_aligned')

aligned_cols

cell_dat <- run.flowsom(cell_dat, aligned_cols, meta.k = 15)

# Dimension reduction -----------------------------------------------------
cell_sub <- do.subsample(cell_dat, targets = 1e4)
cell_sub <- run.umap(cell_sub, aligned_cols)

# Plot results ------------------------------------------------------------
make.colour.plot(cell_sub,
                 "UMAP_X",
                 "UMAP_Y",
                 col.axis = "FlowSOM_metacluster",
                 col.type = "factor",
                 add.label = TRUE,
                 path = "plots")

make.multi.plot(cell_sub, 
                "UMAP_X", 
                "UMAP_Y", 
                plot.by = aligned_cols,
                path = "plots")

make.multi.plot(cell_sub, 
                "UMAP_X", 
                "UMAP_Y", 
                plot.by = "FlowSOM_metacluster", 
                divide.by = "Group", 
                col.type = "factor",
                dot.size = 0.5,
                path = "plots")

# Expression heatmap -----------------------------------------------------
medians <- do.aggregate(cell_dat, aligned_cols, by = "FlowSOM_metacluster")

make.pheatmap(medians, 
              sample.col = "FlowSOM_metacluster", 
              plot.cols = aligned_cols,
              path = "plots")

# Annotate clusters -------------------------------------------------------
annots <- list("T cells"             = c(1, 4),
               "CD45-"               = c(6, 8, 10),
               "Immature B cells"    = 3,
               "Mature B cells"      = 2,
               "Monocytes"           = c(12, 14),
               "Mature neutrophils"  = c(11, 13),
               "Immature neutrophils"= 9,
               "Other/Unknown"   = c(5, 7, 15)
               )

annots <- do.list.switch(annots)
names(annots) <- c("Values", "Population")
# setorderv(annots, 'Values')
annots

cell_dat <- do.add.cols(cell_dat, 
                        base.col = "FlowSOM_metacluster", 
                        add.dat = annots, 
                        add.by = "Values")
cell_dat

cell_sub <- do.add.cols(cell_sub, 
                        base.col = "FlowSOM_metacluster", 
                        add.dat = annots, 
                        add.by = "Values")
cell_sub

make.colour.plot(cell_sub,
                 "UMAP_X",
                 "UMAP_Y",
                 col.axis = "Population",
                 col.type = "factor",
                 add.label = TRUE,
                 path = "plots")

# Plot annotated heatmap --------------------------------------------------
medians <- do.aggregate(cell_dat, transformed_cols, by = "Population")

make.pheatmap(medians, "Population", transformed_cols, path = "plots")

# Export fcs files --------------------------------------------------------
write.files(cell_dat,
            file.prefix = "output/",
            divide.by = "Sample",
            write.csv = FALSE,
            write.fcs = TRUE)

# Summary data ------------------------------------------------------------
simple_names <- c("CD3e", "CD16-32", "Ly6G", "CD45", "CD48", "CD11b", "B220", "Ly6C")
names(cell_dat)[22:29] <- simple_names

sum_dat <- create.sumtable(dat = cell_dat,
                           sample.col = "Sample",
                           pop.col = "Population",
                           use.cols = simple_names,
                           annot.cols = c("Group", "Batch"))
sum_dat

names(sum_dat)

# Plot heatmap comparisons ------------------------------------------------
## Z-score calculation
plot_cols <- names(sum_dat)[4:75]
sum_dat_z <- do.zscore(sum_dat, plot_cols, replace = TRUE)

make.pheatmap(sum_dat_z,
              sample.col = "Sample",
              plot.cols = plot_cols,
              is.fold = TRUE,
              annot.cols = c("Group", "Batch"),
              dendrograms = 'column',
              path = "plots")


# Principal component analysis of samples ---------------------------------
run.pca(sum_dat, 
        use.cols = plot_cols,
        plot.ind.group = TRUE,
        group.ind = "Group",
        path = "PCA"
)

# Statistical analysis ----------------------------------------------------
stats_tab <- create.stats(sum_dat, 
                          use.cols = plot_cols,
                          sample.col = "Sample",
                          group.col = "Group",
                          comparisons = list(c("Mock", "Virus")),
                          corrections = "fdr")

FC <- unlist(stats_tab[1, c(-1, -2)])
p_values <- unlist(stats_tab[2, c(-1, -2)])

stats_df <- data.frame(FC = FC, p_values = p_values, vars = plot_cols)

make.volcano.plot(p_values,
                  FC, 
                  vars = as.character(1:72), 
                  title = "Mock vs Virus",
                  xlim = c(-4, 4),
                  path = "plots")

fold_filter <- abs(stats_df$FC) >= 0.26
p_filter <- stats_df$p_value < 0.05

stats_filtered <- stats_df[fold_filter & p_filter, ]$vars
stats_filtered <- stats_filtered[!is.na(stats_filtered)]

make.pheatmap(sum_dat_z,
              sample.col = "Sample",
              plot.cols = stats_filtered,
              is.fold = TRUE,
              annot.cols = c("Group", "Batch"),
              dendrograms = 'column',
              path = "plots")

make.pheatmap(sum_dat_z,
              sample.col = "Sample",
              plot.cols = stats_filtered,
              is.fold = TRUE,
              annot.cols = c("Group", "Batch"),
              dendrograms = 'both',
              path = "plots")
