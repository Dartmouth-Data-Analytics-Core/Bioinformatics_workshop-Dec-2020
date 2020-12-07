# Basic programming with R

R is a free, open source programming language and statistical software environment that is used extensively in bioinformatics. In addition to the large collection of functions included in base distribution, there are an enormous number of packages designed to extend R's functionality for specific applications.

R is also a very powerful way to create high quality graphics, using both functionality in base R as well as graphics specific packages such as [ggplot2](https://ggplot2.tidyverse.org/). These packages provide a high level of user control, meaning almost all plotting features can be controlled. Importantly, numerous R packages provide functionality for generating bioinformatics specific visualizations.

Visit the *R-Project* homepage [here](https://www.r-project.org/)

![](../figures/r-logo.png)

> **Note:** This is not a comprehensive introduction to R-programming, but rather a review of the basics to help get you started. In addition to the materials provided to you before the workshop, there are some links to more comprhensive tutorials for R in the 'Useful_links.md' in the parent directory of the workshop repository.

## RStudio

RStudio is an IDE ()

![](../figures/r-studio-logo.png)






## Basic data structures in R

Functional programming in R is achieved by assigning values to objects. Objects are specific data structures that take on a particular form defined by that objects *class*.

### Vectors

The most basic object class in R are vectors, and can generally only hold one type of data (a property referred to as being *atomic*).

In R, five basic object classes exist:
- numeric - real numbers (e.g. 1.342)
- integer - whole numbers (e.g. 1.0)
- character - strings of characters (e.g. letters, words, sentences)
- logical - `TRUE` or `FALSE`
- complex - numbers with 'imaginary' parts (not commonly used)

Vectors can be created using the `c()` function, which will combine its arguments together, and the assignment operator `<-` which tells R you want to assign that vector to a specific value.
```r
# numeric
x <- c(1.63, 2.25, 3.83, 4.99)

# integer
x <- c(1, 2, 3, 4)

# characters
x <- c("a", "b", "c", "d")

# logical
x <- c(TRUE, FALSE, TRUE, TRUE)
```

Each object class has specific *attributes*, which we can extract using the appropriate accessor functions. For example, the class of an object is itself an attribute that can be obtained using the `class()` function:
```r
class(x)
```

Another important attribute is length. For example, if we want to know how many elements are in a character string, we can use the `length()` function.
```r
length(x)
```

Vectors can be combined or nested to create a single vector, or evaluated against each other:
```r
# combine a vector and a nested vector
x <- c(1, 2, 3, 4, c(1, 2, 3, 4))

# multiple two integer vector
y <- c(2, 2, 2, 2)
x * y
```

Even though vectors are atomic, they can be coerced from one class to another using functions written to modify their attributes. e.g.
```r
x <- c(1, 2, 3, 4)
as.character(x)

x <- c(TRUE, FALSE, TRUE, TRUE)
as.numeric(x)

x <- c(1.63, 2.25, 3.83, 4.99)
as.integer(x)
```

Elements within vectors can be subset or indexed based on their position in that vector. Individual elements can also be assigned names, which can also be used to perform indexing.
```r
# define a chacter string
x <- c("a", "b", "c", "d")

# get elements 1 and 3
x[c(1,3)]

# get elements 1 to 3
x[c(1:3)]

# define a numeric vector
x <- c(1.63, 2.25, 3.83, 4.99)

# assign it names
names(x) <- c("gene 1", "gene 2", "gene 3", "gene 4")

# index for specific element
x[]"gene 1"]
```

Vectors can contain missing values, defined by `NA` and `NaN`. These elements can be identified with the functions `is.na()` or `is.nan()`:
```r
x <- c(1.63, NA, 3.83, 4.99)
x

x.na <- is.na(x)
x.na

# what object class is returned  
class(x.na)
```

### Factors







Factors are very important when you start performing any statistical analyses in R. 



### Lists

subsetting lists

### Matrices




### Data frames

colnames and rownames
subsetting matricies and dataframes




logical operators (!, &, |, ==)


object classes and data structures/ classes introduced by packages
s3 and as4..?


## Import and export data






## Loops and functions


for()




## Basic data visulization



there are a large number of specific packages designed for visualization. We won't cover these here since they are covered extensively elsewhere, but some useful visualixation packages to be aware of include:
- [ggplot2]()
- [ggpubr]()
- [ploty]()

## Save data in R session

write.table()
write.csv()

save()
load()


saveRDS()
readRDS()
