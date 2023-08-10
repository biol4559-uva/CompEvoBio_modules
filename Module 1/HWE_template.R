### libraries
  .libPaths(c("/project/biol4559-aob2x/biol4559-R-packages/", .libPaths()))
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

### generate data
  ### this will work
    test_fun(x=5, y=2)

  ### so will this
    test_fun(x=c(1:10), y=c(-1:-10))

  ### this helps you keep track of input and output
    out <- test_fun2(x=c(1:10), y=c(-1:-10))

### plot data
  ggplot(data=out, aes(x=x, y=y, color=z)) + geom_point()
