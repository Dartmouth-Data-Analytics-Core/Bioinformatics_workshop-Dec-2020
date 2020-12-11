#04-intro to R
##############################################
print(1, 2, 3)

?print()

# where are you on your local machine
getwd()

# set working directory to the data folder of Day-2 in the github repo you downloaded - notice that the path needs to be in double quotes
setwd("your_path/Bioinformatics_workshop/Day-2/data/")

############################################## Vectors
# numeric
x <- c(1.63, 2.25, 3.83, 4.99)

# integer
x <- c(1, 2, 3, 4)

# characters
x <- c("a", "b", "c", "d")

# logical
x <- c(TRUE, FALSE, TRUE, TRUE)

class(x)

length(x)

# combine a vector and a nested vector
x <- c(1, 2, 3, 4, c(1, 2, 3, 4))
x

# multiply two integer vectors
y <- c(2, 2, 2, 2)
x * y

x <- c(1, 2, 3, 4)
as.character(x)

x <- c(TRUE, FALSE, TRUE, TRUE)
as.numeric(x)

x <- c(1.63, 2.25, 3.83, 4.99)
as.integer(x)

# define a chacter string
x <- c("a", "b", "c", "d")

# get elements 1 and 3 
x[c(1,3)]

# get elements 1 to 3 using the ':' operator
x[c(1:3)]

# define a numeric vector
x <- c(1.63, 2.25, 3.83, 4.99)

# assign it names
names(x) <- c("gene 1", "gene 2", "gene 3", "gene 4")

# index for specific element
x["gene 1"]

x <- c(1.63, NA, 3.83, 4.99)
x

x.na <- is.na(x)
x.na

# what object class is returned  
class(x.na)


############################################## Operators
# addition
1 + 2

# multiplication
2 * 3

# exponentiation
2^4

x <- c(1, 2, 3, 4)

# which elements are less than 3
x < 3

# which elements are less than or equal to 3
x <= 3

# define a character string
x <- c("a", "b", "c", "d", "a")

# which elements are equal to a
x == "a"

x <- c(1, 2, 3, 4)

# which elements are NOT equal to 4
x != 4

# which could also be achieved with
!x == 4

# which elements are less than 2 or equal to 4
x < 2 | x ==4

# which elements are less than 2 AND equal to 4
x < 2 & x ==4

x <- c(1, 2, 3, 4)

# subset x for values less than 3
x_sub <- x[x < 3]
x_sub

# define a character string
x <- c("a", "b", "c", "d", "a")

# subset x for elements equal to a
x[x == "a"]

############################################## Factors
# make a character vector with only male or female as entries
x <- c("female", "female", "male", "female", "male")

# use factor() constructor function to generate the factor
x <- factor(x, levels = c("female", "male"))

# confirm the class and check the levels
class(x)
levels(x)

# use table() to count up all the instances of each level
table(x)

############################################## Lists
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

x <- list(c(1.63, 2.25, 3.83, 4.99),
          c(TRUE, FALSE, TRUE, TRUE),
          c("a", "b", "c", "d"))

# use the structure function str()to examine the composition of the list
str(x)


############################################## Matrices
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

rownames(mat1) <- c("gene_1", "gene_2", "gene_3")
colnames(mat1) <- c("subject_1", "subject_2", "subject_3")

############################################## Data Frames
df <- data.frame(subject_id = c("s1", "s2", "s3", "s4"),
                 age = c(45, 83, 38, 23),
                 gender = c("female", "female", "male", "female"),
                 status = c("case", "case", "control", "control"))

str(df)

df <- data.frame(subject_id = c("s1", "s2", "s3", "s4"),
                 age = c(45, 83, 38, 23),
                 gender = c("female", "female", "male", "female"),
                 status = c("case", "case", "control", "control"),
                 stringsAsFactors=FALSE)

str(df)

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

# obtain a logical indicating which subjects are female
df$gender == "female"

# use logical to subset the data frame for only female subjects (rows)
df2 <- df[df$gender == "female", ]

# check dimensions of the new data frame
dim(df2)

# use the LOGICAL NOT operator ! to obtain only male subjects  
df[!df$gender == "female", ]

# this could obviously also be achieved with..
df[df$gender == "male", ]

############################################## Functions
# define the function that takes a single argument and returns a single argument
myfun <- function(x){
  y <- x + 1
  return(y)
}

# call the function
myfun(x = 10)

# assign the output to a new variable
var1 <- myfun(x = 10)
var1

myfun2 <- function(num1, num2){
  num3 <- num1 + num2
  return(num3)
}

# call the function
myfun2(num1 = 10, num2 = 11)


############################################## Loops
x <- c(1,2,3,4,5,6,7,8,9)

# define the loop, using [i] to define the elements of x used in each iteration
for(i in 1:length(x)){
  print(x[i] * 10)
}

x <- c(1,2,3,4,5,6,7,8,9)

# create variable to save results to
y <- c()

# define and run the loop
for(i in 1:length(x)){
  y[i] <- x[i] * 10
}

return(y)

############################################## Data Visualization
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

# generate a simple histogram for x
hist(x, col = "indianred")

# the breaks argument can be used to change how the intervals are broken up
hist(x, col = "indianred", breaks=10)
hist(x, col = "indianred", breaks=50)
hist(x, col = "indianred", breaks=100)

############################################## Data import/export
# using read.table
counts <- read.table(file = "all_counts.txt",
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)
### Note 1: header accepts logical value indicating if the first row are column names (default FALSE)
### Note 2: we use stringsAsFactors

# check class, dimensions and structure
class(counts); dim(counts); str(counts)

# using read.delim
counts <- read.table(file = "all_counts.txt", stringsAsFactors=FALSE)

# check class, dimensions and structure
class(counts); dim(counts); str(counts)

# subset counts for first 5 columns and first 2000 genes
counts_sub <- counts[1:2000, 1:5]

# write to tab delimited file using write.table
write.table(counts_sub, file = "all_counts_sub.txt", sep = "\t")

write.csv(counts_sub, file = "all_counts_sub.csv")

############################################## Saving R objects
# create some R objects
x <- c(1.63, 2.25, 3.83, 4.99)
y <- c(TRUE, FALSE, TRUE, TRUE)
z <- c("a", "b", "c", "d")

# save all 3 objects to one file 
save(x, y, z, file = "my_r_objects.rdata")

load(file = "my_r_objects.rdata")

# save a single object to a specific file path
saveRDS(x, file = "my_r_object.rds")

# use the assignment operator to assign object to any variable you want
x <- readRDS(file = "my_r_object.rds")

# I changed my mind, I want the object to be assigned to variable `a` in my new env
a <- readRDS(file = "my_r_object.rds")
