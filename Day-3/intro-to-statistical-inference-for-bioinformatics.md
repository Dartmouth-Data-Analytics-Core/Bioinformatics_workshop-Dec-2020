
# Introduction to Statistical Learning & Inference in bioinformatics 

As we discussed on day 1, bioinformatics draws on knowledge from multiple disciplines. To effectively solve most bioinformatic problems, knowledge from each of these disciplines must be applied to address the overall problem. Developing a **working knowledge of statistics** is critical for anyone hoping to solve bioinformatic problems, particularly in genomics. 

<p align="center">
  <img src="../figures/bioinfo-venn.png" height="400" width="400"/>
</p>

Adapted from Md Arman Hossen on [Medium](https://medium.com/datadriveninvestor/i-have-designed-my-own-bioinformatics-degree-260b24767d87). 

Statistics involves the analysis of numerical data that we wish to use for making inferences on a larger population. Generally, ***statistical learning*** refers to the models and procedures with which we analyze these numerical datasets. 

These models and procedures can be used to either make **predictions** (e.g. healthy vs. diseased tissue based on gene expression profiles) or for ***statistical inference***, where tools from statistical learning are applied to help us understand the relationship between a set of variables (e.g. which genes are associated with diseased tissue). 

While a comprehensive introduction to *statistical learning and inference* is well beyond the scope of this workshop (and generally takes years of specific training and experience), we will introduce some fundamental concepts required for approriate statsitcal analysis of genomics data sets that are relevant to wet- or dry-lab scientists alike. 

**What we *will* cover:**  
- What is statistical learning? 
- Basics of supervised & unsupervised data analysis methods
- Fundamentals of statistical inference & hypothesis testing 
- *P*-values & multiple testing correction 

**What we *will not* cover:**  
- Math & probability theory underlying statistical learning procedures
- Model selection procedures & how to assess model fit 
- A comprehensive introduction to the methods used for statistical learning and inference in bioinformatics

> **Important note:** Builing a more complete understanding of the statistical procedures commonly used in bioinformatics, such that you are able to confidently implement, interpret, and troubleshoot these procedures, requires a strong working knowledge of relevant math and probability theory. Such training is best sought out through formal instruction, and is usually not included in applied bioinformatics courses. While developing an understanding of the fundamental concepts in statistical learning and inference covered here will allow you to begin leveraging more complex statistical analyses in your own research, and act as a solid fondation upon which to further your training in this domain, it is also important to recognize when more specialist expertise is needed to address your analytical questions. 

--- 

## Statistical learning 

As described above, *statistical learning* describes the models and procedures we use in order to understand data. Generally, the methods we use can be categorized into *supervised* and *unsupervised* methods. 

**Supervised methods** describe approaches used when a set of observations we have made (e.g. gene expression levels) are all associated with a response variable (e.g. diseased or healthy). The observations, also called predictors, are generally referred to as the *independent variable (X)*, while the response is the *dependant variable (Y)*. 

By fitting statistical models to the data (the *learning* part), we aim to learn about the relationship between the *independent* and *dependent* variables (the *inference* part). For example, we may apply a form of linear model to gene expression data from heathly cases and diseased controls in order to address questions like:  
- *Which genes are associated with disease?*
- *How much does a given gene contribute to the disease phenotype?*

**Examples of supervised methods:**   
- Linear regression
- Generalized linear regression models 
- Descision trees 
- Support vector machines 

**Unsupervised methods** 

There are times when observations are not associated with a predictor (the *dependent variable*) and we simply wish to explore the relationships that exist in our data in a way that is *not supervised* by any such variable. This is often true in *exploratory data analysis* when we want to explore how samples are related to each other without assiging samples to a specific group for modeling purposes. For example, we may wish to confirm samples in an RNA-seq experiment fo not cluster by batch. 

Alternatively, we may not have a dependent variable that can be used to model our observations. For example, in analysis of single cell sequencing data, we are often interested in studying subpopulations of cells that come from the same sample, and therefore require some way assessing similarities and differences between the cells so that wse can identify potential subpopulations of interest.  

**Examples of unsupervised methods:**    
- Dimensionality reduction (e.g. PCA, NMF, t-SNE, UMAP)
- Clustering methods (e.g. hierachical, K-means)
- Hidden markov modeling 

Below, we provide more speicifc introductions to both supervised and unsupercised learnings, using basic linear modeling as an example for supervised approaches, while exploring PCA and hierachical clustering for unsupervised analysis. 

> A more comprehensive introduction to statistical learning can be found in the book: [An Introduction to Statistical Learning](http://faculty.marshall.usc.edu/gareth-james/ISL/). 






## Statistical inference 


statsitical inference refers to: 
estimation
hypothesis testing 


mention MLE 

bayseian 


### What is *hypothesis testing*?







### What is a *P*-value?





### The multiple testing problem






### Methods for multiple testing correction













### Supervised learning - Linear modeling 

Simple linear models, or linear regression, is used pervasively in bioinformatics and genomics for statistical inference. Linear models are relatively simple, flexible, and interpretable, meaning they make excellent tools for statistical inference and scale well to thousands of observations, which is critical for common genomics datasets. Example applications of linear models include:  
- RNA-seq (differential expression)
- ChIP-seq (differential binding)
- ATAC-seq (differential accessibility)
- Microarray analysis (e.g. DNA methylation)
- Variant identification (WES/WGS/RNA-seq)
- Genome-wide association studies (GWAS)

Understanding the basics of linear modeling is central to being able to perform these types of analyses in a statistical programming environment such as R. 

Given their importance and pervasive use in bioinformatics and genomics, we will introduce the fundamental concepts of linear models, and how you can fit these models in R. 

> **Note:** Linear modeling is the topic of entire independent courses and again requires knowledge of appropriate mathematics and proprability to understand completely. Thus, this should be considered an introduction rather than a standalone resurce. 

In a standard linear model, we assume that some *response* variable (*Y*) can be represented as a linear combination of a set of *predictors* (*X*, independent variables). In building a linear model, we estimate a set of *coefficients* that explain how the *predictors* are related to the *response*. We can use these *coefficients* for statistical inference to better understand which predictors are associated with our response, or for applying the model to new data where we wish to predict the response variable using only a set of predictors. 

Before reviewing out the statistical notation for a simple linear model, it can be useful to first consider the main components: 
`*response* = *predictor(s)* + *error*`

- **The *response*** is the dependent variable we wish to model based on some set of predictors  

- **The *predictor(s)*** is the independent variable, or variables, that we wish to model as a linear combination of the response (this is usually what we are interested in for statiatical inference and hypothesis testing)  

- **The *error*** component represents the information not explained by the model, and exists since we know that no set of predictors will perfectly explain the response variable. These are often referred to as the *residuals* of the model.  

Using the statistical notation for a simple linear regression model:

Y = &beta;<sub>0</sub> +  &beta;<sub>i</sub> X<sub>i</sub> + &epsilon;

- Y is a continuous response variable that we assume in normally distributed
- &beta;<sub>i</sub> are the coefficients to be estimated (&beta;<sub>i</sub>-value)
- X<sub>i</sub> are the predictors 
- &beta;<sub>0</sub> refers to the model intercept
- &epsilon; refers to the error term (residuals) and are assumed to follow a normal distribution

There can be any (reasonable) number of predictors (X) in a model, and can be either *continuous* (e.g. age) or categorical (e.g. treatment group, batch). 

Each predictor is associated with a coefficient that describes the relationship of that predictor to the response variable. In the context of linear regression, the coefficient is also referred to as the *slope*. 

In R, the basic syntax for this model is: `lm(response ~ predictor)`. Lets simulate some data that we can use to illustrate the theory described above and fit out first linear model in R. 

```r 
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
``` 

![]("../figures/lm_example-0.png")

The coefficient for the independent/predictor variable, Hba1c, describes its relation to the response variable, expression of gene X. Here, the coefficient in telling us that *for every 1 unit increase in gene expression measured, Hba1c levels increase by ~0.96 units*. 

This is basic *statistical inference*, as we have used this procedure to model the relationship between two variables, and *infer* something about how those variables are related. 

To help us better understand the model, we can plot the regression line on our scatterplot. 
```r
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
```

![]("../figures/lm_example.png")

The regression line (shown in black) illustrates the clear linear relationship between exprerssion of gene X and Hba1c levels. 

The residuals (blue lines) describe how far away each observation (the gene expression values) are from the predicted values from the linear model. All observations are close to the regression line, suggesting the model is a good fit for the data. 

**However**, by virtue of this being a statistical model, all coefficients are estimated with some level of uncertainty. If the model is a poor fit for the data, there will be a high uncertainty in the coefficient.

One way to evaluate how much meaning we should attribute to the coefficient, is to calculate a *P*-value for it through hypothesis testing, which we will explore below. 

> **Note:** Although standard models for modeling gene expression data would include expression values as the response variable, these models usually take on a more complicated form (see note on *Generalized linear models* below), however we have set up a more simple model for teaching purposes. 

#### Hypothesis testing with linear models 

In order to test how much certainty we have for a particular coefficient from a linear model, we estimate a quantity called **the standard error (SE)**. Without discussing the underlying stastics that define it, the SE is essentially a *measure of certainty around the coefficient*, and is dependent on the variance of the residuals (&epsilon;). 

Importantly, the SE can be used to perform **hypothesis testing** to determine if the coefficient is statistically significant. In this case, we can test the null hypothesis that the coefficient is equal to zero, using the following equation to calculate the *t-score*: 

*t-score* = &beta;<sub>i</sub> - 0 / SE(&beta;<sub>i</sub>)

The *t-score* can then be used to calculate a *P*-value, as described in the hypothesis testing section. In R, the `summary()` function will test all model coefficients against the null hypothesis: 
```r
sum_lm1 <- summary(lm1)
sum_lm1

# get the coefficients table 
coef(sum_lm1)

# get the coefficients themselves 
coef(sum_lm1)[,1]

# get the P-value for the hba1c coefficient 
coef(sum_lm1)[2,4]
```

The *P*-value is very small, so we can reject the null, and conclude that Hba1c levels are associated with expression of gene X, and interpret the coefficient as a meaningful quantity. 

If the *P*-value does not pass the *a priori*significance threshold for your analysis, the coefficient should be ignored as that predcitor is **not associated** with the response variable. 

You can always confirm by looking at the slope in a simple linear model. To demonstrate this, explore the example below for Gene Y and its relation to Hba1c levels. 
```r
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

```


![]("../figures/lm_example-2.png")


The flatter slope of the regression line, and larger values of the residuals, suggests there is no useful relationship between Hba1c levels and expression of gene Y, which is supported by the large *P*-value returned by the model. 


#### Simple Linear modeling with categorical variables 

In genomics, we commonly have categorial predictor variables, in contrast to the continuous variable (Hba1c) from our example above. Example of categorial variable include: 
- Wild-type vs knockout
- Vehicle vs treatment 
- Control vs diseased

Importantly, linear models are capable of incorportaing categorical variables as predictors. Lets consider another example, where we have gene expression levels for gene X measured in 20 healthy tissues, and 20 diseased tissues, and we wish to use a linear model to explore the relation between gene expression and disease status. 

```r
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
abline(lm1, lty=2)
```

![]("../figures/lm_example-3.png")

Looking at the model output, the *P*-value is very small, therefore we can conclude that there is an association between expression of gene X and disease status in this sample. 

Again, the coefficient tells us about the relationship between the predictor and the response. The coefficient for the predictor `subject_group` tells us that for each unit increase in this variable, there is an increase of 11.2 expression units for gene X. 

Since a 'unit increase' in `subject_group` simply means controls vs diseased subjects, we can interpret this as the difference in expression between controls and cases. This is analgous to how we would calculate a fold-change value in an RNA-seq analysis. 

---
**Multiple regression**

We could have simply addressed the above analysis using a more simple statistical test such as a *t-test*. However, we commonly want to include additional variables in our linear models, and approaches such as the t-test cannot handle this scenario. 

For example, we might want to control for factors that might confound gene expression differences between the control and diseased groups. For example, we could control for age and sex of the subjects, or perhaps the batch the samples were collected in. 

In this scenario, we can use linear models to control for the additional variables by adding them into the statsitical model e.g.
```r
lm(dat3$exp_geneX ~ dat3$subject_group + dat3$age + dat3$age + dat3$age)
```

This approach is referred to as **multiple regression**. If you will be doing any sort of complex bioinformatics data analysis involving linear models, I strongly encourage you to use this primer as a starting point to learn more about multiple regression and more complex linear modeling scenarios.  

---

#### Generalized linear models

Commonly in bioinformatics, we find ourselves needing a model that can assume a different statistical distribution than the normal. The rationale behind this is that many bioinformatics and genomic data types exhibit specific distributions that are different from the normal. 

We fit the data using a **generalized linear model (GLM)**. GLM's are a family of statistical models that generalize standard linear regression in two ways:  
- use of probability distributions other than the normal distribution 
- the use of a *link-function* that links the expression values in the linear model to the experimental groups, in a way that these other distributions can be used. 

![]("../figures/neg-binom.png")

For example, bulk RNA-seq data typically exhibit a distribution referred to as the *negative-binomial* and therefore require a GLM of the *negative-binomial family* in order to appropriately model RNA-seq counts and test them for differential expression. 

While GLMs are beyond the scope of this workshop, and we simply do not have the time to cover them in this short course, we do cover the fundamentals of how GLMs are used in the context of RNA-seq data analysis in our RNA workshop! 



### Unsupervised learning - Dimension reduction & clustering 

As mentioned above, two commonly used types of unsupervised learning are *dimension reduction* and *clustering-based* methods. Both encompass a number of distinct methodologies that have various strengths and weaknes


# use CB book to help get an example function to use 
```r
pheatmap()
```


As mentioned above, two commonly used types of unsupervised learning are *dimension reduction* and *clustering-based* methods. Both encompass a number of distinct methodologies that have various strengths and weaknesses. 

To gain an appreciation for how these methods are used and presented in bioinformatic and genomic data analyses, we will explore the fundamental aspects of *principal components analysis (PCA)* as an example of dimension reduction, and *hierachical clustering* as an example of clustering-based analyses. 

#### Principal components analysis 



t-SNE and UMAP 



```r
prcomp()
```






