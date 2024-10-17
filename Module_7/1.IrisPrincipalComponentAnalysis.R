
####################################################################
########## First PCA section of Iris flower dataset ################
####################################################################

### install some new packages: poolfstat, ggforce, zoo, ggrepel

# Libraries
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggforce) # install.packages("ggforce")
  library(zoo) # install.packages("zoo")
  library(ggrepel) # install.packages("ggrepel")

# Load Iris flower dataset
  data(iris)

# Convert to data.frame object
  iris.dt <- data.frame(iris)

# YOUR TURN: make the data from wide --> long

# Plot data in facets
  ggplot(data=irisL[!data=="Species"],
         aes(x=value,
             fill=Species)) +
    geom_histogram(bins=30) +
    facet_wrap(~data, nrow = 2) +
    theme_bw() +
    labs(x="Value (mm)",
         y="Number of samples") +
    theme(strip.text = element_text(face="bold.italic", size=16),
          title = element_text(size=18, face="bold"),
          legend.text = element_text(size=16, face="bold.italic"),
          legend.title = element_text(face="bold", size=18),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))

# Introduce geom_mark_eclipse
  ggplot(data=iris.dt,
         aes(x=Petal.Length,
             y=Sepal.Length,
             color=Species,
             fill=Species)) +
    geom_point() +
    geom_mark_ellipse() +
    theme_bw() +
    labs(x="Petal length (mm)",
         y="Sepal length (mm)") +
    theme(strip.text = element_text(face="bold.italic", size=16),
          title = element_text(size=18, face="bold"),
          legend.text = element_text(size=16, face="bold.italic"),
          legend.title = element_text(face="bold", size=18),
          axis.text.x = element_text(face="bold", size=18),
          axis.text.y = element_text(face="bold", size=18),
          axis.title.x = element_text(face="bold", size=20),
          axis.title.y = element_text(face="bold", size=20),
          axis.title = element_text(face="bold", size=20))

# Normalize the data using scale() function with centering and scaling
  normalized_data <- scale(iris.dt[,-5], center = TRUE, scale = TRUE)

# Perform PCA on the normalized data
  iris.pca <- prcomp(normalized_data)
  summary(iris.pca)

# look inside the iris.pca object
  str(iris.pca)
  pca.dt <- data.table(PC1=iris.pca$x[,1], PC2=iris.pca$x[,2], species=iris.dt$Species)
  pca.dt

# Your turn: Scatter plot of the first two PCs colored by species

# How much variation does the PCs contribute to?
# Lets make a scree plot to investigate.
summ.pca <- summary(iris.pca)
prop.var <- data.table(var=c(summ.pca$importance[2,]*100),
                       PC=colnames(summ.pca$importance))

# Basic scree plot
  barplot(prop.var$var, names.arg = prop.var$PC,
          xlab = "Principal Component",
          ylab = "Proportion of Variance (%)",
          main = "Proportion of Variance Explained by Principal Components",
          col = "blue")

# Utilize both the data.table, ggplot2, ggforce, & patchwork packages to add color/labels/eclipses

# New ggforce package and eclipses function
# Add: geom_mark_ellipse() to the ggplot object to get those cool looking eclipses around each species!
