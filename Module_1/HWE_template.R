### specify your library path
  .libPaths(c("~/USERNAME", .libPaths())) ### <- you'll need to change USERNAME to your computing ID

### install packages; you only need to this once (ever)
  install.packages("ggplot2", "data.table")

### libraries
  library(ggplot2)
  library(data.table)

### function
  test_fun <- function(x, y) {
    return(x*y)
  }

  test_fun2 <- function(x, y) {
    z <- x*y
    return(data.table(x=x, y=y, z=z))
  }

  test_fun3 <- function(x,y) {
    y1 <- x^2
    y2 <- x^1.5

    return(data.table(x=x, y1=y1, y2=y2))
  }

### generate data
### this function returns a single value
  test_fun(x=5, y=2)

### the same function can return a vector of values.
  test_fun(x=c(1:10), y=c(-1:-10))

### does this work? What does the error message mean?
  test_fun(x=c(1:10), y=c(1:3))

### the output of this function helps you keep track of input and output by returning a data.table
  out <- test_fun2(x=c(1:10), y=c(-1:-10))
  ggplot(data=out, aes(x=x, y=y, color=z)) + geom_point()

### there are a few basic data types in R:
  a <- "alan" ### this is vector of length 1
  b <- c("alan", "robert") ### this is a vector of length 2. We need to use the "c" function to concatenate the two single values
  c <- data.table(name=c("alan", "robert"), role=c("teacher", "TA"), number=c(2,1)) ### this is a data.table/data.frame. It is a 2 dimensional array, with two columns and two rows. Each column can take a different type of data, like a character string, a number, etc
  
### This function performs two operations on the input value and returns a WIDE table. We can plot y1 and y2 on the same graph in a few differnet ways
    out2 <- test_fun3(x=c(1:100))

    ### method 1: manually code the two y-values.
    ggplot(data=out2) +
      geom_line(aes(x=x, y=y1), color="blue") +
      geom_line(aes(x=x, y=y2), color="orange")

    ### convert wide to long and use groupings in ggplot
    out2_long <- melt(out2, id.vars="x", measure.vars=c("y1", "y2"))

    ggplot(data=out2_long) +
      geom_line(aes(x=x, y=value, group=variable, color=variable))

### saving plots is easy
  p1 <- ggplot(data=out2_long) +
    geom_line(aes(x=x, y=value, group=variable, color=variable))

  ggsave(p1, file="~/test_plot1.png")
