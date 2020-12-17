#05-Intro to stats
########################################

# generate an expty vector to store P-values in
p.value <- c()

# generate 2 random variables (r.v.s) 1000 times and run a t.test on each pair of r.v.s
# the r.v.s will have the same mean and distribution
for(i in 1:1000){
  # simulate random variables
  x <- rnorm(n = 20, mean = 0, sd = 1)
  y <- rnorm(n = 20, mean = 0, sd = 1)
  # run t.test and extract p-value
  p.value[i] <- t.test(x, y)$p.value
}

# count number of P-values less than our alpha
table(p.value < 0.05)

# order vector of P-values
p.value <- p.value[order(p.value)]

# visualize P-value magnitude
plot(-log10(p.value), las = 1,
     col = "cornflowerblue",
     main = "-log 10 P-values",
     xlab = "ordered P-values")
#### Note: taking -log10 of small number gives big number!

# add a horizontal line
abline(h=-log10(0.05), col = "red", lty = 2)

# add some useful labels
text(600, 1.5, "Small P-values higher up")
text(600, 1.1, "Large P-values lower down")
########################################

# visualize P-value magnitude
plot(-log10(p.value), las = 1,
     col = "cornflowerblue",
     main = "-log 10 P-values",
     xlab = "ordered P-values", ylim = c(0,5))

# add a horizontal line
abline(h=-log10(0.05), col = "red", lty = 2)
text(600, 4.5, "Bonferroni")

# add a horizontal line
abline(h=-log10(0.05/1000), col = "red", lty = 2)
text(600, 1.5, "Original threshold")
########################################

# bonferroni correction
p.adj.bonf <- p.adjust(p.value, method = "FWER")
p.adj.bonf

# check if any are signifciant
table(p.adj.bonf < 0.05)
########################################

#BiocManager::install("qvalue")
# library(qvalue)
qvalue(p.value)$qvalues

# check how many are sig.
table(p.adj.fdr < 0.05)
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
boxplot(dat3$exp_geneX ~ dat3$subject_group ,
     ylab = "Expression (Gene X)",
     xlab = "Subject group",
     main = "Gene X exp. vs Hba1c",
     col = c("indianred", "cornflowerblue"), pch = 16, las = 1)


# run the linear model and evaluate
lm_2 <- lm(dat3$exp_geneX ~ dat3$subject_group)
summary(lm_2)
