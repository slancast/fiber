#/usr/bin/R
#This will be to test the factor analysis and the linear descriminant analysis

#First up the factor analysis
#From https://www.promptcloud.com/blog/exploratory-factor-analysis-in-r/

data <- read.csv("/Users/SLancaster/Downloads/EFA.csv",header=TRUE)
head(data)

library('psych')
#We’ll be using `Psych` package’s `fa.parallel` function to execute parallel analysis. Here we specify the data frame and factor method (`minres` in our case). Run the following to find acceptable number of factors and generate the `scree plot`:
parallel <- fa.parallel(data, fm = 'minres', fa = 'fa')

#Now that we’ve arrived at probable number number of factors, let’s start off with 3 as the number of factors. In order to perform factor analysis, we’ll use `psych` package’s  `fa()function. 
threefactor <- fa(data,nfactors = 3,rotate = "oblimin",fm="minres")
print(threefactor)
fa.diagram(threefactor)

#Now for LDA
#From 
#Already problems in the frst line, the website calls the data wdbc.csv but it's rather called
#wdbc.data. Makes it very difficult to figure out what is going on
library(MASS)

data(iris)

r <- lda(formula = Species ~ ., 
         data = iris, 
         prior = c(1,1,1)/3)

r$prior
r$counts
r$means
r$scaling
r$svd

prop = r$svd^2/sum(r$svd^2)

r2 <- lda(formula = Species ~ ., 
          data = iris, 
          prior = c(1,1,1)/3,
          CV = TRUE)

head(r2$class)
head(r2$posterior, 3)

train <- sample(1:150, 75)

r3 <- lda(Species ~ ., # training model
          iris, 
          prior = c(1,1,1)/3, 
          subset = train)

plda = predict(object = r, # predictions
               newdata = iris[-train, ])

head(plda$class)
head(plda$posterior, 3)
head(plda$x, 3)



