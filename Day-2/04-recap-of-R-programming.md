# Basic programming with R

R is a free, open source programming language and statistical software environment that is used extensively in bioinformatics. In addition to the large collection of functions included in base distribution, there are an enormous number of packages designed to extend R's functionality for specific applications.

R is also a very powerful way to create high quality graphics, using both functionality in base R as well as graphics specific packages such as [ggplot2](https://ggplot2.tidyverse.org/). These packages provide a high level of user control, meaning almost all plotting features can be controlled. Importantly, numerous R packages provide functionality for generating bioinformatics specific visualizations.

Visit the *R-Project* homepage [here](https://www.r-project.org/)

<img src="../figures/r-logo.png" height="110" width="150"/>

> **Note:** This is not a comprehensive introduction to R-programming, but rather a review of the basics to help get you started. In addition to the materials provided to you before the workshop, there are some links to more comprhensive tutorials for R in the 'Useful_links.md' in the parent directory of the workshop repository.

## RStudio

RStudio is an IDE ()

<img src="../figures/r-studio-logo.png" height="130" width="300"/>






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

# get elements 1 and 3 using the ':' operator
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

### Operators

We introduced two *operators* in the examples above, the assignment operator `<-` and the sequence operator `:`. Operators are essentially symbols that tell R how you would like to relate the *operands* on either side of the symbol. in R, operators can be broadly categorized into *assignment*, *arithmetic*, *relational*, and *logical*.

The major assignment operators are `<-` and `=` which both tell R to assign a vector to a some value. It is generally safer to use `<-` (in my opinion).

Below are the basic *arithmetic*, *relational*, and *logical* operators that are useful to know.

**Arithmetic operators**

Operator | Effect
----- | -----
+ | addition
-	| subtraction
*	| multiplication
/	| division
^ |	exponentiation

Some example usage of arithmetic operators:
```r
# addition
1 + 2

# multiplication
2 * 3

# exponentiation
2^4
```

**Relational operators**

Operator | Effect
----- | -----
<	| less than
>	| greater than
<= |	less than or equal to
>=	| greater than or equal to
== |	equal to
== |	Note equal to

Some example usage of relational operators:
```r
x <- c(1, 2, 3, 4)

# which elements are less than 3
x < 3

# which elements are less than or equal to 3
x <= 3

# define a character string
x <- c("a", "b", "c", "d", "a")

# which elements are equal to a
x == "a"
```

**Logical operators**

Operator | Effect
----- | -----
!	| NOT
&	| AND
\|	| OR

Some example usage of logical operators:
```r
x <- c(1, 2, 3, 4)

# which elements are NOT equal to 4
x != 4

# which could also be achieved with
!x == 4

# which elements are less than 2 and equal to 4
x < 2 | x ==4
```

Relational and logical operators can be used to subset a vector based on the values returned by the operator, and the brackets, as we did above for specific elements.
```r
x <- c(1, 2, 3, 4)

# subset x for values less than 3
x_sub <- x[x < 3]


# define a character string
x <- c("a", "b", "c", "d", "a")

# subset x for elements equal to a
x[x == "a"]
```


### Factors

Factors are a special instance of vectors where only predefined values, called *levels* can be included in the vector. Such vectors are useful when you know that elements of a vector should take on one of those predfined values.

Categorical data is often stored in vectors, making them a very important object class when you start doing any statistical modeling in R. For example, you might store subject gender for all the subjects in your study as a factor, with the levels *male* and *female*.

```r
# make a character vector with only male or female as entries
x <- c("female", "female", "male", "female", "male")

# use factor() constructor function to generate the factor
x <- factor(x, levels = c("female", "male"))

# confirm the class and check the levels
class(x)
levels(x)

# use table() to count up all the instances of each level
table(x)
```

### Lists

Sometimes, it may be desirable to store multiple vectors, or even vectors of different object classes, into the same R overall object. Lists are a special object class that permits objects with these attributes, making them distinct from atomic vectors.

In the same way that vectors and factors are constructed using `c()` and `factors()` respectively, lists are created using the `lists()` constructor function.

```r
x <- list(c(1.63, 2.25, 3.83, 4.99),
          c(2.43, 8.31, 3.12, 7.25),
          c(1.29, 3.23, 3.48, 0.23))

# the structure function str() can be useful to examine the composition of a list
str(x)

# confirm the length
length(x)

# lists can be subset using brackets
### subset for first element of list
x[[1]]
### subset for first element of first vector in list
x[[1]][1]

# lists elements can be given names using a character vector equal to list length
names(x) <- c("gene_1", "gene_2", "gene_3")

# names can also be used to subset a list
x[["gene_1"]]

# subsetting can also be achieved with the $ subsetting operator
x$gene_1
```

In our example list, all three vectors stored in the list are numeric, however as mentioned above, lists can store vectors of different classes.


```r
x <- list(c(1.63, 2.25, 3.83, 4.99),
          c(TRUE, FALSE, TRUE, TRUE),
          c("a", "b", "c", "d"))

# the structure function str() can be useful to examine the composition of a list
str(x)
```

### Matrices

By extending the attributes of a vector to give them *dimensions*, i.e. the number of rows and columns we want the vector to be organized into, we can create *matrices*, a data structure that efficiently stores tabular data of a specific, single object class.

```r
mat <- matrix(c(1.63, 2.25, 3.83, 4.99),
              c(2.43, 8.31, 3.12, 7.25),
              c(1.29, 3.23, 3.48, 0.23),
              nrow=3, ncol=4)
# check the structure and dimensions with dim()
str(mat)
dim(mat)

# specific elements can be obtained through subsetting
### row 1
mat1[1,]
### column 2
mat1[,2]
### element 2 of row 3
mat[3,2]

# check class of the object and one row
class(mat)
class(mat[1,])
```

Since matrices have dimensions, `names()` cannot be used as we did for vectors. Instead, `names()` is generalized into `rownames()` and `colnames()`.

```r
rownames(mat1) <- c("gene_1", "gene_2", "gene_3")
colnames(mat1) <- c("subject_1", "subject_2", "subject_3")
```

Matrices are a very important object class for mathematical and statistical applications in R, so it is certainly worth exploring more complex matrix operations if you will be doing any more complex statistical analysis in R.

### Data frames

Data frames are very efficient ways of storing tabular data in R. Like matrices, data frames have dimensionality and are organized into rows and columns, however data frames can store vectors of different object classes.

Often you will construct a data frame by reading in a dataset from a file. While we will cover reading data into R below, we will construct a data frame using the `data.frame()` constructor function in R for this example.

```r
df <- data.frame(subject_id = c("s1", "s2", "s3", "s4"),
                 age = c(45, 83, 38, 23),
                 gender = c("female", "female", "male", "female"),
                 status = c("case", "case", "control", "control"))

str(df)
```

Note that the default behavior of `data.frame()` is to convert character strings to factors. If you want to prevent this behavior, you can set the `StringsAsFactors()` argument as `FALSE`.

```r
df <- data.frame(subject_id = c("s1", "s2", "s3", "s4"),
                 age = c(45, 83, 38, 23),
                 gender = c("female", "female", "male", "female"),
                 status = c("case", "case", "control", "control"),
                 stringsAsFactors=FALSE)

str(df)
```

Data frames can be subset in similar ways to matrices using brackets or the `$` subsetting operator. Columns/variables can also be added using the `$` operator.
 ```r
# get first row
df[1,]

# get first column
df[,1]

# get gender variable/column
df[, c("gender")]

# # get gender and status
df[, c("gender", "status")]

# get the gender variable with $
df$gender

# add a column for smoking status
df$smoking_status <- c("former", "none", "heavy", "none")
```

Relational (e.g. `==`) and logical operators (e.g. `!`) can be used to interrogate specific variables in a dataframe. The resulting logical can also be used to subset the data frame.
```r
# obtain a logical indicating which subjects are female
df$gender == "female"

# use logical to subset the data frame for only female subjects (rows)
df2 <- df[df$gender == "female", ]

# check dimensions of the new data frame
dim(df2)

# use the LOGICAL NOT operator ! to obtain only male subjects  
df[!df$gender == "female", ]

# this could obviously also be achieved with..
df[!df$gender == "male", ]
```

### Beyond the basic object classes

As we discussed, one of the major advantages of R is its enormous user base that are continuously developing and releasing packages. Implementing the additional functionality of these packages often requires **more complex data object classes to be created**, which are generally related in some way to one or more of the basic data structures in R that we have discussed.

The general approach used to create these novel classes is referred to as *object-orientated programming (OOP)*. Although we will not go into any detail on OOP in this workshop, it is useful to know that several OOP methods are used to create classes in R.

The [S4 OOP approach](http://adv-r.had.co.nz/S4.html) is perhaps the most relevant to this workshop, as it is heavily used by packages from the [Bioconductor project](http://bioconductor.org/), which we will be using on Day 3. S4 OOP provides a flexible approach to creating objects with multiple *slots*, each with defined names and object classes.

An example of a package that has used an S4 OOP approach to create objects of novel classes is the Bioconductor package [*GenomicRanges*](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), which provides representation and query of genomic region data in R through object classes such as `GRanges`.

To efficiently represent genomic regions data, `GRanges` class objects, at least 3 key slots (each with their own associated class for that vector) will be needed:
- a chromosome slot: with class `character` (e.g. chr1)
- a start coordinate slot: with class `integer` (e.g. 338833)
- an end coordinate slot: with class `integer` (e.g. 338987)

Creating object classes in this way is desirable as it allows *accessor functions* to be created, which allow very simple interaction with the objects of this class. For example, simply using the `chr()` accessor function to easily extract all chromosome identities of the genomic regions in my object.

An excellent tutorial describing S4 classes and their use in the [Bioconductor project](http://bioconductor.org/) is available [here](https://bioconductor.org/help/course-materials/2017/Zurich/S4-classes-and-methods.html). While this is not a topic you need understand in great detail, it is worth understanding the basic principles.


## Functions

Beyond the functions implemented in base R and packages that you install, R allows you to create user defined functions, which can perform any functionality that you can define.

Defining your own functions can be useful when you want to perform a specific set of tasks repeatedly on some input(s) and return a defined output. Furthermore, once defined functions become part of your global environment and are therefore preserved for future use, minimizing the need for repetitive code.

Functions are created using using `function()` with the assignment operator `<-`. The arguments you use in the `function()` command define the variables that those arguments will be assigned to when you call the function. The last line of the function defines what output is returned.

Let's define a basic function as an example.
```r
# define the function
myfun <- function(x){
  y <- x + 1
  return(y)
}

# call the function
myfun(x = 10)

# assign the output to a new variable
var1 <- myfun(x = 10)
var1
```

Functions can have as many arguments as you specify. The names of these arguments are only assigned as variables within the function, however it is good practice to avoid using arguments with the same name as variables already existing in your global environment.

For example, if I already have a variable named `x` in my environment, I should avoid using x to define the name of the arguments to my function.
```r
myfun2 <- function(num1, num2){
  num3 <- num1 + num2
  return(num3)
}

# call the function
myfun2(num1 = 10, num2 = 11)
```

## Loops

Loops are used to iterate over a piece of code multiple times, and can therefore be used to achieve specific tasks. The most often used type of loop in R is the `for()` loop, which will evaluate the contents of the loop for all of the values provided to the `for()` function.

For example:
```r
x <- c(1,2,3,4,5,6,7,8,9)

# define the loop, using [i] to define the elements of x used in each iteration
for(i in 1:length(x)){
  print(x[i] * 10)
}
```

We may wish to save the output of each iteration of the loop to a new variable, which can be achieved using the following approach:
```r
x <- c(1,2,3,4,5,6,7,8,9)

# create variable to save results to
y <- c()

# define and run the loop
for(i in 1:length(x)){
  y[i] <- x[i] * 10
}

return(y)
```

## Basic data visualization

R is a very powerful tool for visualization and provides a large amount of user control over the plotting attributes. Basic visualization in R is achieved using the `plot()` function.

```r
# generate a set of random numbers to plot
x <- rnorm(1000, mean = 10, sd = 2)
y <- rnorm(1000, mean = 20, sd = 1)

# plot x against y to produce a scatter plot
plot(x, y)

# add labels, title, and color
plot(x, y,
     main = "X vs. Y",
     xlab = "X values",
     ylab = "Y values",
     col = "red")
```

R can also be easily used to generate histograms:
```r
# generate a simple histogram for x
hist(x, col = "indianred")

# the breaks argument can be used to change how the intervals are broken up
hist(x, col = "indianred", breaks=10)
hist(x, col = "indianred", breaks=50)
hist(x, col = "indianred", breaks=100)
```

### Visualization specific packages

There are a large number of packages designed for specifically for visualization and are very useful in bioinformatic analyses. We won't cover these here since they are covered extensively elsewhere, but some useful visualization packages to be aware of include:
- [ggplot2]()
- [ggpubr]()
- [ploty]()

Importantly, visualization implemented in these packages form the basis for some bioinformatics specific data visualization packages that we will explore later in the workshop.

## Import and export data








## Save data in R session

write.table()
write.csv()

save()
load()


saveRDS()
readRDS()
