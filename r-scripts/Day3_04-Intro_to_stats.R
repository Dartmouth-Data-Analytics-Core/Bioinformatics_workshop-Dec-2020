## 03-Intro to Stats

########################################
# read in the example data
dat <- read.csv("lm-example-data.csv", stringsAsFactors=FALSE)

# explore it quickly
head(dat)
str(dat)

# plot
plot(dat$gene_exp ~ dat$hba1c,
     ylab = "Expression (Gene X)",
     xlab = "Hba1c score",
     main = "Gene X exp. vs Hba1c",
     col = "indianred", pch = 16, las = 1)

# fit a linear model with gene expression as the response
lm1 <- lm(dat$gene_exp ~ dat$hba1c)
lm1

########################################
# generate plot again
plot(dat$gene_exp ~ dat$hba1c,
     ylab = "Expression (Gene X)",
     xlab = "Hba1c score",
     main = "Gene X exp. vs Hba1c",
     col = "indianred", pch = 16, las = 1)

# add the model on the scatterplot
abline(lm1, lty=2)

# calculate the predicted gene expression values using the model
pre <- predict(lm1)

# plot the difference between the predicted and the true values
segments(dat$hba1c, dat$gene_exp, dat$hba1c, pre,
         col="cornflowerblue")
#### Note: These are the residuals!

########################################
sum_lm1 <- summary(lm1)
sum_lm1

# get the coefficients table
coef(sum_lm1)

# get the coefficients themselves
coef(sum_lm1)[,1]

# get the P-value for the hba1c coefficient
coef(sum_lm1)[2,4]

########################################
# read in the example data
dat2 <- read.csv("lm-example-data-geneY.csv", stringsAsFactors=FALSE)

# plot
plot(dat2$gene_exp ~ dat2$hba1c,
     ylab = "Expression (Gene Y)",
     xlab = "Hba1c score",
     main = "Gene Y exp. vs Hba1c",
     col = "indianred", pch = 16, las = 1)

# fit a linear model with gene expression as the response
lm1 <- lm(dat2$gene_exp ~ dat2$hba1c)
summary(lm1)

# add the model on the scatterplot
abline(lm1, lty=2)

# plot the difference between the predicted and the true values
segments(dat2$hba1c, dat2$gene_exp, dat2$hba1c, pre, col="cornflowerblue")

########################################
# read in the example data
dat3 <- read.csv("lm-example-3.csv", stringsAsFactors=FALSE, row.names = 1)

# quickly explore it
head(dat3)
table(dat3$subject_group)
# Note: Controls are coded as 0, cases are coded as 1

# visualize the data
plot(dat3$subject_group, dat3$exp_geneX,
     ylab = "Expression (Gene X)",
     xlab = "Subject group",
     main = "Gene X exp. vs Hba1c",
     col = "indianred", pch = 16, las = 1)


# run the linear model and evaluate
lm_2 <- lm(dat3$exp_geneX ~ dat3$subject_group)
summary(lm_2)

# add regression line to the plot
abline(lm_2, lty=2)
## this isn't plotting a line on the plot -need to change this to lm_2 instead of lm1

######################################## 
# read in data
fpkm <- read.table("fpkm_sub.txt", stringsAsFactors = FALSE, header=TRUE)
meta <- read.delim("metadata_sub.tsv", stringsAsFactors = FALSE, header = TRUE)

# have a quick look at top of files
head(fpkm[,1:6])
head(meta)

# log transform fpkm counts
log_fpkm <- log2(fpkm+1)

# calculate variance across samples for each gene
vars <- apply(log_fpkm, 1, var)
# order variances based on magnitude of variance
vars <- rev(vars[order(vars)])

# plot variance for genes accross samples
plot(vars, las = 1, main="Sample gene expression variance", xlab = "Gene", ylab = "Variance")
abline(v=5000, col="red")

######################################## PCA
# perform PCA and order by variance
vars_sub <- vars[1:5000]

# perform the PCA on the fpkm matrix
pca <- prcomp(t(log_fpkm[names(vars_sub), ]))

# look at the object returned
str(pca)
head(pca$x)

# construct data frame w/ PC loadings and add sample labels
pca_df <- as.data.frame(pca$x)
pca_df$tissue <- as.factor(meta$Biosample.term.name)
pca_df$sample_ids <- meta$File.accession

# extract percent variance explained by each PC
percentVar <- pca$sdev^2/sum(pca$sdev^2)

# add colors for plotting to df
cols <- grDevices::rainbow(length(levels(pca_df$tissue)))

# create an empty variable in pca_df to be filled with colors from cols
pca_df$col <- NA

# loop over tissue types in pca_df and assign colors for plotting
for(i in 1:length(levels(pca_df$tissue))){
  ind1 <- which(pca_df$tissue == levels(pca_df$tissue)[i])
  pca_df$col[ind1] <- cols[i]
}

# plot PC1 vs PC2
plot(pca_df[, 1], pca_df[, 2],
     xlab = paste0("PC1 (", (round(percentVar[1], digits=3)*100), "% variance)"),
     ylab = paste0("PC2 (", (round(percentVar[2], digits=3)*100), "% variance)"),
     main = paste0("PC1 vs PC2 genome-wide RNA-seq tissue profiles"),
     pch = 16, cex = 1.35, cex.lab = 1.3, cex.axis = 1.15, las = 1,
     panel.first = grid(),
     col = pca_df$col)

# add a legend to the plot
legend(1.5, 105, levels(pca_df$tissue), pch = 16, col = cols, cex = 0.9)

######################################## Hierarchical Clustering
# subset mnetadata to 5 tissues of interest
meta_ord <- meta[meta$Biosample.term.name=="forebrain" |
                   meta$Biosample.term.name=="heart" |
                   meta$Biosample.term.name=="limb" |
                   meta$Biosample.term.name=="liver" |
                   meta$Biosample.term.name=="intestine", ]

# subset FPKM matrix to contain the same subset of samples
log_fpkm_sub <- log_fpkm[, c(colnames(log_fpkm) %in%  meta_ord$File.accession)]

# calculate variance of each gene across samples for new subset of data
vars <- apply(log_fpkm_sub, 1, var)

# order variances based on magnitude of variance
vars <- rev(vars[order(vars)])

# plot variance for genes across samples
plot(vars, las = 1, main="Sample gene expression variance",
     xlab = "Gene", ylab = "Variance")
# add vertical line
abline(v=1000, col="red")

########################################
# subset var to only top 2000 genes with most variance
vars_sub <- vars[1:2000]

# subset the fpkm matrix to these genes
mat <- log_fpkm_sub[names(vars_sub), ]

# order the samples in the same order they are in in the metadata file
mat <- mat[, c(match(meta_ord$File.accession, colnames(mat)))]

# scale the fpkm matrix by row
mat_scaled = t(apply(mat, 1, scale))

# set column names for this matrix (they were removed during transposition)
colnames(mat_scaled) <- colnames(mat)

########################################
# load the pheatmap package
library(pheatmap)

# create data frame to annotate heatmap with
annotation_col = data.frame(Tissue = meta_ord$Biosample.term.name)
rownames(annotation_col) = meta_ord$File.accession

# use pheatmap() to perform clustering on scaled data matrix
pheatmap(mat_scaled,
         show_rownames=FALSE, show_colnames=FALSE,
         annotation_col = annotation_col,
         cluster_cols = TRUE,
         clustering_method = "average",
         clustering_distance_cols = "correlation")
